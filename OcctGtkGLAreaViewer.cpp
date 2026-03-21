#include "OcctGtkGLAreaViewer.h"
#include "OcctGlTools.h"
#include "OcctGtkTools.h"
#include "OcctInputBridge.h"
#include <Message.hxx>
#include <OpenGl_Context.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <OpenGl_FrameBuffer.hxx>

#if !defined(__APPLE__) && !defined(_WIN32) && defined(__has_include)
  #if __has_include(<Xw_DisplayConnection.hxx>)
    #include <Xw_DisplayConnection.hxx>
    #define USE_XW_DISPLAY
  #endif
#endif
#ifndef USE_XW_DISPLAY
typedef Aspect_DisplayConnection Xw_DisplayConnection;
#endif

#ifdef _WIN32
#else
  #include <gdk/x11/gdkx.h>
#endif

#include <AIS_ListOfInteractive.hxx>
#include <AIS_Shape.hxx>
#include <AIS_ViewCube.hxx>

OcctGtkGLAreaViewer::OcctGtkGLAreaViewer()
{
  Handle(Aspect_DisplayConnection) aDisp = new Xw_DisplayConnection();
  Handle(OpenGl_GraphicDriver) aDriver = new OpenGl_GraphicDriver(aDisp, false);

  aDriver->ChangeOptions().buffersNoSwap = true;
  aDriver->ChangeOptions().buffersOpaqueAlpha = true;
  aDriver->ChangeOptions().useSystemBuffer = false;
  aDriver->ChangeOptions().contextCompatible = false;

  myViewer = new V3d_Viewer(aDriver);
  myViewer->SetDefaultBackgroundColor(Quantity_NOC_WHITE);
  myViewer->SetDefaultLights();
  myViewer->SetLightOn();
//  myViewer->ActivateGrid(Aspect_GT_Rectangular, Aspect_GDM_Lines);

// =========================================================================
  // NEW: Configure TOPMOST layer to ignore Z-buffer (depth)
// =========================================================================
  {
      // Get current settings for the topmost layer
      Graphic3d_ZLayerSettings aSettings = myViewer->ZLayerSettings(Graphic3d_ZLayerId_Topmost);

      // Disable Depth Test (to draw over everything)
      aSettings.SetEnableDepthTest(false);

      // Disable Depth Write (so it doesn't affect subsequent drawing)
      aSettings.SetEnableDepthWrite(false);

      // Clear depth before drawing this layer
      aSettings.SetClearDepth(true);

      // Apply settings back
      myViewer->SetZLayerSettings(Graphic3d_ZLayerId_Topmost, aSettings);
  }
// =========================================================================

  myContext = new AIS_InteractiveContext(myViewer);

  myViewCube = new AIS_ViewCube();
  myViewCube->SetViewAnimation(myViewAnimation);
  myViewCube->SetFixedAnimationLoop(false);
  myViewCube->SetAutoStartAnimation(true);

  myView = myViewer->CreateView();
  myView->SetImmediateUpdate(false);
  myView->ChangeRenderingParams().ToShowStats = true;
  myView->ChangeRenderingParams().StatsPosition = new Graphic3d_TransformPers(
      Graphic3d_TMF_2d, Aspect_TOTP_RIGHT_LOWER, Graphic3d_Vec2i(5, 5));
  myView->ChangeRenderingParams().CollectedStats = Graphic3d_RenderingParams::PerfCounters_FrameRate;
  // NOLINTNEXTLINE
 /// myView->ChangeRenderingParams().CollectedStats = (Graphic3d_RenderingParams::PerfCounters)(
 ///   Graphic3d_RenderingParams::PerfCounters_FrameRate | Graphic3d_RenderingParams::PerfCounters_Triangles);

#ifdef HAVE_GLES2
  set_use_es(true);
#endif

  signal_realize()  .connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaRealized));
  signal_unrealize().connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaReleased), false);
  signal_render()   .connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaRender), false);

  // Initialize Bridge
  m_input_bridge = std::make_unique<OcctInputBridge>(
      this,
      myView,
      myContext,
      this,
      [this]() { this->queue_draw(); },

      // TEXT CALLBACK: Send signal to the window
      [this](const std::string& msg) {
          if (!msg.empty()) {
              signal_status_message.emit(msg);
          }
      }
  );

  AIS_MouseGestureMap& aMouseGestures = ChangeMouseGestureMap();
  aMouseGestures.Clear();
  aMouseGestures.Bind(Aspect_VKeyMouse_MiddleButton, AIS_MouseGesture_RotateOrbit);
  aMouseGestures.Bind(Aspect_VKeyMouse_MiddleButton | Aspect_VKeyFlags_CTRL, AIS_MouseGesture_RotateOrbit);
  aMouseGestures.Bind(Aspect_VKeyMouse_RightButton, AIS_MouseGesture_Pan);
  aMouseGestures.Bind(Aspect_VKeyMouse_LeftButton, AIS_MouseGesture_SelectRectangle);

  AIS_MouseSelectionSchemeMap& aMouseSelScheme = ChangeMouseSelectionSchemes();
  aMouseSelScheme.Clear();
  aMouseSelScheme.Bind(Aspect_VKeyMouse_LeftButton, AIS_SelectionScheme_Replace);
  aMouseSelScheme.Bind(Aspect_VKeyMouse_LeftButton | Aspect_VKeyFlags_SHIFT, AIS_SelectionScheme_XOR);
}

OcctGtkGLAreaViewer::~OcctGtkGLAreaViewer() {}

bool OcctGtkGLAreaViewer::onModifiersChanged(Gdk::ModifierType theType) { return false; }

void OcctGtkGLAreaViewer::updateModifiers() {
  Aspect_VKeyFlags aFlags = Aspect_VKeyFlags_NONE;
  if (AIS_ViewController::Keys().IsKeyDown(Aspect_VKey_Shift)) aFlags |= Aspect_VKeyFlags_SHIFT;
  if (AIS_ViewController::Keys().IsKeyDown(Aspect_VKey_Control)) aFlags |= Aspect_VKeyFlags_CTRL;
  if (AIS_ViewController::Keys().IsKeyDown(Aspect_VKey_Meta)) aFlags |= Aspect_VKeyFlags_META;
  if (AIS_ViewController::Keys().IsKeyDown(Aspect_VKey_Alt)) aFlags |= Aspect_VKeyFlags_ALT;
  myKeyModifiers = aFlags;
}

bool OcctGtkGLAreaViewer::onKeyPressed(guint theKeyVal, guint theKeyCode, Gdk::ModifierType ) {
  const Aspect_VKey aVKey = OcctGtkTools::gtkKey2VKey(theKeyVal, theKeyCode);
  if (aVKey == Aspect_VKey_UNKNOWN) return false;
  const double aTimeStamp = AIS_ViewController::EventTime();
  AIS_ViewController::KeyDown(aVKey, aTimeStamp);
  updateModifiers();
  AIS_ViewController::ProcessInput();
  return true;
}

void OcctGtkGLAreaViewer::onKeyReleased(guint theKeyVal, guint theKeyCode, Gdk::ModifierType ) {
  const Aspect_VKey aVKey = OcctGtkTools::gtkKey2VKey(theKeyVal, theKeyCode);
  if (aVKey == Aspect_VKey_UNKNOWN) return;
  const double aTimeStamp = AIS_ViewController::EventTime();
  AIS_ViewController::KeyUp(aVKey, aTimeStamp);
  updateModifiers();
  AIS_ViewController::ProcessInput();
}

void OcctGtkGLAreaViewer::onMotionMove(double theX, double theY) {
  if (OcctGtkTools::gtkHandleMotionEvent(*this, myView, Graphic3d_Vec2d(theX, theY), myKeyModifiers)) queue_draw();
}

void OcctGtkGLAreaViewer::onMouseButtonPressed(int , double theX, double theY) {}
void OcctGtkGLAreaViewer::onMouseButtonReleased(int , double theX, double theY) {}
bool OcctGtkGLAreaViewer::onMouseScroll(double theDeltaX, double theDeltaY) { return true; }

void OcctGtkGLAreaViewer::handleViewRedraw(const Handle(AIS_InteractiveContext)& theCtx, const Handle(V3d_View)& theView) {
  AIS_ViewController::handleViewRedraw(theCtx, theView);
  if (myToAskNextFrame) queue_draw();
}

void OcctGtkGLAreaViewer::initPixelScaleRatio() {
  AIS_ViewController::SetTouchToleranceScale(myDevicePixelRatio);
  myView->ChangeRenderingParams().Resolution = (unsigned int )(96.0 * myDevicePixelRatio + 0.5);
}

void OcctGtkGLAreaViewer::dumpGlInfo(bool theIsBasic, bool theToPrint) {
  TColStd_IndexedDataMapOfStringString aGlCapsDict;
  myView->DiagnosticInformation(aGlCapsDict, theIsBasic ? Graphic3d_DiagnosticInfo_Basic : Graphic3d_DiagnosticInfo_Complete);
  TCollection_AsciiString anInfo;
  for (TColStd_IndexedDataMapOfStringString::Iterator aValueIter(aGlCapsDict); aValueIter.More(); aValueIter.Next()) {
    if (!aValueIter.Value().IsEmpty()) {
      if (!anInfo.IsEmpty()) anInfo += "\n";
      anInfo += aValueIter.Key() + ": " + aValueIter.Value();
    }
  }
  myGlInfo = anInfo;
  if (theToPrint) Message::SendInfo(anInfo);
}

void OcctGtkGLAreaViewer::onGlAreaRealized() {
  make_current();
  Graphic3d_Vec2i aLogicalSize(get_width(), get_height());
  if (aLogicalSize.x() == 0 || aLogicalSize.y() == 0) return;

  try {
    OCC_CATCH_SIGNALS
    throw_if_error();
    Glib::RefPtr<Gdk::Surface> aGdkSurf = this->get_native()->get_surface();
    myDevicePixelRatio = get_scale_factor();
    initPixelScaleRatio();
    const bool isFirstInit = myView->Window().IsNull();
    Aspect_Drawable aNativeWin = 0;
#if defined(_WIN32)
#else
    aNativeWin = gdk_x11_surface_get_xid(aGdkSurf->gobj());
#endif
    const Graphic3d_Vec2i aViewSize = Graphic3d_Vec2i(Graphic3d_Vec2d(aLogicalSize) * myDevicePixelRatio + Graphic3d_Vec2d(0.5));
    if (!OcctGlTools::InitializeGlWindow(myView, aNativeWin, aViewSize, myDevicePixelRatio)) {
        // error handling
    }
    make_current();
    dumpGlInfo(true, true);
    if (isFirstInit) myContext->Display(myViewCube, 0, 0, false);
  } catch (...) {}
  myContext->SetPixelTolerance(3);
}

void OcctGtkGLAreaViewer::onGlAreaReleased() {
  make_current();
  try {
    throw_if_error();
    Handle(Aspect_DisplayConnection) aDisp;
    if (!myContext.IsNull()) {
      aDisp = myViewer->Driver()->GetDisplayConnection();
      myContext->RemoveAll(false);
      myContext.Nullify();
      myView->Remove();
      myView.Nullify();
      myViewer.Nullify();
    }
    make_current();
    aDisp.Nullify();
  } catch (...) {}
}

bool OcctGtkGLAreaViewer::onGlAreaRender(const Glib::RefPtr<Gdk::GLContext>& theGlCtx) {
  (void )theGlCtx;
  if (myView->Window().IsNull()) {
    onGlAreaRealized();
    make_current();
    gtk_gl_area_attach_buffers(gobj());
  }
  if (myView->Window().IsNull()) return false;

  try {
    throw_if_error();
    if (!OcctGlTools::InitializeGlFbo(myView)) return false;
    const Graphic3d_Vec2i aLogicalSize(get_width(), get_height());
    Graphic3d_Vec2i aViewSize; myView->Window()->Size(aViewSize.x(), aViewSize.y());
    const float aPixelRatio = float(aViewSize.y()) / float(aLogicalSize.y());
    if (myDevicePixelRatio != aPixelRatio) {
      myDevicePixelRatio = aPixelRatio;
      Handle(OcctGlTools::OcctNeutralWindow) aWindow = Handle(OcctGlTools::OcctNeutralWindow)::DownCast(myView->Window());
      aWindow->SetDevicePixelRatio(aPixelRatio);
      initPixelScaleRatio();
      dumpGlInfo(true, false);
    }
    myView->InvalidateImmediate();
    AIS_ViewController::FlushViewEvents(myContext, myView, true);
    return true;
  } catch (...) { return false; }
}

static std::string shapeTypeToString(TopAbs_ShapeEnum type) {
    switch (type) {
        case TopAbs_COMPOUND:  return "COMPOUND";
        case TopAbs_COMPSOLID: return "COMPSOLID";
        case TopAbs_SOLID:     return "SOLID";
        case TopAbs_SHELL:     return "SHELL";
        case TopAbs_FACE:      return "FACE";
        case TopAbs_WIRE:      return "WIRE";
        case TopAbs_EDGE:      return "EDGE";
        case TopAbs_VERTEX:    return "VERTEX";
        case TopAbs_SHAPE:     return "SHAPE";
        default:               return "UNKNOWN";
    }
}

static void processXCAFLabelRecursive(
    const TDF_Label& label,
    const Glib::RefPtr<Gtk::TreeStore>& store,
    const Gtk::TreeNodeChildren& children,
    const Handle(XCAFDoc_ShapeTool)& shapeTool,
    const Gtk::TreeModelColumn<Glib::ustring>& colName,
    const Gtk::TreeModelColumn<Glib::ustring>& colType,
    const Gtk::TreeModelColumn<int>& colId,
    std::map<int, TDF_Label>& labelMap,
    std::map<int, TopLoc_Location>& locMap,
    std::map<int, int>& depthMap,
    std::map<int, TopoDS_Shape>& subShapeMap,
    int& idCounter,
    const TopLoc_Location& parentLoc,
    int currentTreeDepth)
{
    TopLoc_Location myLocalLoc = shapeTool->GetLocation(label);
    TopLoc_Location myGlobalLoc = parentLoc * myLocalLoc;
    TDF_Label definition_label = label;
    if (shapeTool->IsReference(label)) shapeTool->GetReferredShape(label, definition_label);

    Handle(TDataStd_Name) aNodeName;
    Glib::ustring name_str = "Unnamed";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName))
        name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();
    else if (definition_label.FindAttribute(TDataStd_Name::GetID(), aNodeName))
        name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();

    TopoDS_Shape partShape;
    shapeTool->GetShape(definition_label, partShape);
    std::string typeStr = "Unknown";
    if (!partShape.IsNull()) typeStr = shapeTypeToString(partShape.ShapeType());

    auto row = *store->append(children);
    row[colName] = name_str;
    row[colType] = typeStr;

    int currentId = idCounter++;
    row[colId] = currentId;

    labelMap[currentId] = definition_label;
    locMap[currentId]    = myGlobalLoc;
    depthMap[currentId] = currentTreeDepth;

    TDF_ChildIterator child_it(definition_label, Standard_False);
    if (child_it.More()) {
        for (; child_it.More(); child_it.Next()) {
            processXCAFLabelRecursive(child_it.Value(), store, row.children(), shapeTool,
                                      colName, colType, colId,
                                      labelMap, locMap, depthMap, subShapeMap,
                                      idCounter, myGlobalLoc, currentTreeDepth + 1);
        }
    }
    else {
        if (!partShape.IsNull() && partShape.ShapeType() == TopAbs_SOLID) {
            TopExp_Explorer exp(partShape, TopAbs_FACE);
            int faceIdx = 1;
            for (; exp.More(); exp.Next()) {
                const TopoDS_Shape& face = exp.Current();
                auto faceRow = *store->append(row.children());
                faceRow[colName] = Glib::ustring::compose("Face %1", faceIdx++);
                faceRow[colType] = "FACE";
                int faceId = idCounter++;
                faceRow[colId] = faceId;
                subShapeMap[faceId] = face.Moved(myGlobalLoc);
            }
        }
    }
}

void OcctGtkGLAreaViewer::activateSubShapeSelection()
{
    if (myContext.IsNull()) return;
    AIS_ListOfInteractive listOfObjects;
    myContext->DisplayedObjects(listOfObjects);
    for (const auto& obj : listOfObjects) {
        if (obj->IsKind(STANDARD_TYPE(AIS_ViewCube))) continue;
        myContext->Deactivate(obj);
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_FACE));
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_EDGE));
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_VERTEX));
    }
}

void OcctGtkGLAreaViewer::selectShapeById(int id)
{
    myContext->ClearSelected(false);
    if (!m_selected_part_ais.IsNull()) {
        myContext->Remove(m_selected_part_ais, false);
        m_selected_part_ais.Nullify();
    }
    bool foundInScene = false;
    if (m_label_map.find(id) != m_label_map.end()) {
        TDF_Label targetLabel = m_label_map[id];
        AIS_ListOfInteractive listOfObjects;
        myContext->DisplayedObjects(listOfObjects);
        for (const auto& obj : listOfObjects) {
            Handle(XCAFPrs_AISObject) xcafObj = Handle(XCAFPrs_AISObject)::DownCast(obj);
            if (!xcafObj.IsNull()) {
                if (xcafObj->GetLabel().IsEqual(targetLabel)) {
                    myContext->AddOrRemoveSelected(obj, false);
                    foundInScene = true;
                    break;
                }
            }
        }
    }
    if (!foundInScene) {
        TopoDS_Shape shapeToSelect;
        if (m_subshape_map.find(id) != m_subshape_map.end()) {
            shapeToSelect = m_subshape_map[id];
        }
        else if (m_label_map.find(id) != m_label_map.end() && m_loc_map.find(id) != m_loc_map.end()) {
             Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
             TopoDS_Shape rawShape;
             ShapeTool->GetShape(m_label_map[id], rawShape);
             if (!rawShape.IsNull()) {
                 shapeToSelect = rawShape.Moved(m_loc_map[id]);
             }
        }
        if (!shapeToSelect.IsNull()) {
            m_selected_part_ais = new AIS_Shape(shapeToSelect);
            myContext->Display(m_selected_part_ais, false);
            myContext->AddOrRemoveSelected(m_selected_part_ais, false);
        }
    }
    myContext->UpdateCurrentViewer();
    queue_draw();
}

void OcctGtkGLAreaViewer::fitToShapeById(int id)
{
    TopoDS_Shape shapeToFit;
    if (m_label_map.count(id) && m_loc_map.count(id)) {
        Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
        TopoDS_Shape rawShape;
        ShapeTool->GetShape(m_label_map[id], rawShape);
        if (!rawShape.IsNull()) {
            shapeToFit = rawShape.Moved(m_loc_map[id]);
        }
    }
    else if (m_subshape_map.count(id)) {
        shapeToFit = m_subshape_map[id];
    }
    if (shapeToFit.IsNull()) return;
    Bnd_Box bbox;
    BRepBndLib::Add(shapeToFit, bbox);
    if (!bbox.IsVoid()) {
        myView->FitAll(bbox, 0.01, false);
        myView->Redraw();
    }
}

int OcctGtkGLAreaViewer::getSelectedShapeId()
{
    myContext->InitSelected();
    if (!myContext->MoreSelected()) return -1;
    Handle(SelectMgr_EntityOwner) owner = myContext->SelectedOwner();
    Handle(StdSelect_BRepOwner) brepOwner = Handle(StdSelect_BRepOwner)::DownCast(owner);
    if (brepOwner.IsNull()) return -1;
    TopoDS_Shape selectedShape = brepOwner->Shape();
    for (auto const& [id, faceShape] : m_subshape_map) {
        if (selectedShape.IsSame(faceShape)) return id;
    }
    if (selectedShape.IsNull()) return -1;
    Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    int bestId = -1;
    int maxDepth = -1;
    for (auto const& [id, label] : m_label_map) {
        if (m_loc_map.find(id) == m_loc_map.end()) continue;
        TopLoc_Location globalLoc = m_loc_map[id];
        TopoDS_Shape partShape;
        ShapeTool->GetShape(label, partShape);
        if (partShape.IsNull()) continue;
        partShape.Move(globalLoc);
        bool isMatch = false;
        if (partShape.IsSame(selectedShape)) {
            isMatch = true;
        } else {
            TopAbs_ShapeEnum selType = selectedShape.ShapeType();
            if (partShape.ShapeType() <= selType) {
                TopTools_IndexedMapOfShape map;
                TopExp::MapShapes(partShape, selType, map);
                if (map.Contains(selectedShape)) isMatch = true;
            }
        }
        if (isMatch) {
            int currentDepth = 0;
            if (m_depth_map.find(id) != m_depth_map.end()) currentDepth = m_depth_map[id];
            if (currentDepth > maxDepth) {
                maxDepth = currentDepth;
                bestId = id;
            }
        }
    }
    return bestId;
}

bool OcctGtkGLAreaViewer::loadStepFile(const Glib::ustring& filePath,
                                Glib::RefPtr<Gtk::TreeStore>& treeModel,
                                const Gtk::TreeModelColumn<Glib::ustring>& colName,
                                const Gtk::TreeModelColumn<Glib::ustring>& colType,
                                const Gtk::TreeModelColumn<int>& colId)
{
    m_is_wireframe = false;
    TCollection_AsciiString path(filePath.c_str());
    Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
    bool isFirstLoad = m_doc.IsNull();
    if (isFirstLoad) {
        Handle(TDocStd_Document) newDoc;
        app->NewDocument("MDTV-XCAF", newDoc);
        if (newDoc.IsNull()) app->NewDocument("XCAF", newDoc);
        m_doc = newDoc;
        myContext->RemoveAll(true);
        if (!myViewCube.IsNull()) myContext->Display(myViewCube, 0, 0, false);
        if (!m_ais_box.IsNull()) m_ais_box.Nullify();
        m_label_map.clear(); m_loc_map.clear(); m_depth_map.clear(); m_subshape_map.clear();
        m_next_id = 1;
        if(treeModel) treeModel->clear();
    }
    Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    TDF_LabelSequence rootsBefore;
    ShapeTool->GetFreeShapes(rootsBefore);
    Standard_Integer countBefore = rootsBefore.Length();
    STEPCAFControl_Reader reader;
    reader.SetColorMode(true); reader.SetNameMode(true); reader.SetLayerMode(true);
    if (reader.ReadFile(path.ToCString()) != IFSelect_RetDone) return false;
    if (!reader.Transfer(m_doc)) return false;
    TDF_LabelSequence rootsAfter;
    ShapeTool->GetFreeShapes(rootsAfter);
    for (Standard_Integer i = countBefore + 1; i <= rootsAfter.Length(); i++) {
        TDF_Label label = rootsAfter.Value(i);
        Handle(XCAFPrs_AISObject) xcafPrs = new XCAFPrs_AISObject(label);
        Handle(Prs3d_Drawer) drawer = xcafPrs->Attributes();
        Handle(Prs3d_ShadingAspect) aspect = drawer->ShadingAspect();
        if (aspect.IsNull()) { aspect = new Prs3d_ShadingAspect(); drawer->SetShadingAspect(aspect); }
        Graphic3d_MaterialAspect mat(Graphic3d_NOM_SHINY_PLASTIC);
        mat.SetSpecularColor(Quantity_NOC_WHITE); mat.SetShininess(0.9);
        aspect->SetMaterial(mat);
        if (xcafPrs->HasColor()) { Quantity_Color color; xcafPrs->Color(color); aspect->SetColor(color); }
        else { aspect->SetColor(Quantity_NOC_GRAY70); }
        myContext->Display(xcafPrs, AIS_Shaded, 0, false);
        myContext->Activate(xcafPrs, AIS_Shape::SelectionMode(TopAbs_EDGE));
        myContext->Activate(xcafPrs, AIS_Shape::SelectionMode(TopAbs_VERTEX));
        myContext->Activate(xcafPrs, AIS_Shape::SelectionMode(TopAbs_FACE));
        if(treeModel) {
            TopLoc_Location startLoc;
            processXCAFLabelRecursive(label, treeModel, treeModel->children(), ShapeTool,
                                      colName, colType, colId,
                                      m_label_map, m_loc_map, m_depth_map, m_subshape_map,
                                      m_next_id, startLoc, 0);
        }
    }
    activateSubShapeSelection();
    myContext->UpdateCurrentViewer();
    myView->FitAll();
    return true;
}