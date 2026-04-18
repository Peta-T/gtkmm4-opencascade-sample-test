#include "OcctGtkGLAreaViewer.h"
#include "OcctGlTools.h"
#include "OcctGtkTools.h"
#include "OcctInputBridge.h"
#include <Message.hxx>
#include <OpenGl_Context.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <OpenGl_FrameBuffer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pln.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Cylinder.hxx>
#include <Graphic3d_AspectFillArea3d.hxx>
#include <Quantity_Color.hxx>
#include <gdk/gdkkeysyms.h>
#include <iostream>
#include <gp_Lin.hxx>
#include <cstdio>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <TDF_Tool.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Prs3d_LineAspect.hxx>
#include <gp_Vec.hxx>
#include <glibmm/main.h>
#include <Geom_Axis2Placement.hxx>
#include <AIS_Trihedron.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <AIS_ColoredShape.hxx>
#include <TopoDS_Iterator.hxx>
#include <cstdlib>
#include <XCAFPrs.hxx>
#include <XCAFPrs_Style.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <filesystem>
#include <Message_ProgressIndicator.hxx>
#include <functional>

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

// Helper class for tracking OpenCASCADE operation progress
class GtkStepProgress : public Message_ProgressIndicator {
private:
    std::function<void(int)> m_callback;
    int m_lastPercent = -1;
public:
    GtkStepProgress(std::function<void(int)> cb) : m_callback(cb) {}

protected:
    virtual void Show(const Message_ProgressScope& theScope, const Standard_Boolean isForce) override {
        // GetPosition() returns a value from 0.0 to 1.0
        int percent = static_cast<int>(GetPosition() * 100.0);

        // Update UI only when percentage changes to avoid slowing down computation
        if (percent != m_lastPercent || isForce) {
            m_callback(percent);
            m_lastPercent = percent;
        }
    }
};

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

  // Configure TOPMOST layer to ignore Z-buffer (depth) for overlays
  {
      Graphic3d_ZLayerSettings aSettings = myViewer->ZLayerSettings(Graphic3d_ZLayerId_Topmost);
      aSettings.SetEnableDepthTest(false);   // Draw over everything
      aSettings.SetEnableDepthWrite(false);  // Don't affect subsequent drawing
      aSettings.SetClearDepth(true);         // Clear depth before drawing this layer
      myViewer->SetZLayerSettings(Graphic3d_ZLayerId_Topmost, aSettings);
  }

  myContext = new AIS_InteractiveContext(myViewer);
  Handle(Prs3d_Drawer) ctxDrawer = myContext->DefaultDrawer();

  // Deflection relative to part size
  ctxDrawer->SetTypeOfDeflection(Aspect_TOD_RELATIVE);
  // Linear deflection (default is 0.001). Smaller means more precise.
  ctxDrawer->SetDeviationCoefficient(0.001);
  // Angular deflection for curved surfaces (e.g., 12 degrees).
  ctxDrawer->SetDeviationAngle(12.0 * M_PI / 180.0);

  // Dynamic highlight color (Orange)
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_Dynamic)->SetColor(Quantity_NOC_ORANGE);
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_LocalDynamic)->SetColor(Quantity_NOC_ORANGE);

  // Selection color
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_Selected)->SetColor(Quantity_NOC_ORANGE);
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_LocalSelected)->SetColor(Quantity_NOC_ORANGE);

  myViewCube = new AIS_ViewCube();
  myViewCube->SetSize(50.0);
  myViewCube->SetFontHeight(10.0);
  myViewCube->SetTransparency(0.8);

  myViewCube->SetViewAnimation(myViewAnimation);
  myViewCube->SetFixedAnimationLoop(false);
  myViewCube->SetAutoStartAnimation(true);

  Handle(Prs3d_DatumAspect) datumAspect = myViewCube->Attributes()->DatumAspect();
  if (datumAspect.IsNull()) {
        datumAspect = new Prs3d_DatumAspect();
        myViewCube->Attributes()->SetDatumAspect(datumAspect);
  }

  // Set colors for axes labels
  datumAspect->TextAspect(Prs3d_DatumParts_XAxis)->SetColor(Quantity_NOC_RED);
  datumAspect->TextAspect(Prs3d_DatumParts_YAxis)->SetColor(Quantity_NOC_GREEN);
  datumAspect->TextAspect(Prs3d_DatumParts_ZAxis)->SetColor(Quantity_NOC_BLUE);

  myView = myViewer->CreateView();
  myView->SetImmediateUpdate(false);
  myView->ChangeRenderingParams().ToShowStats = true;
  myView->ChangeRenderingParams().StatsPosition = new Graphic3d_TransformPers(
      Graphic3d_TMF_2d, Aspect_TOTP_RIGHT_LOWER, Graphic3d_Vec2i(5, 5));
  myView->ChangeRenderingParams().CollectedStats = Graphic3d_RenderingParams::PerfCounters_FrameRate;

#ifdef HAVE_GLES2
  set_use_es(true);
#endif

  signal_realize()  .connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaRealized));
  signal_unrealize().connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaReleased), false);
  signal_render()   .connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onGlAreaRender), false);

  // Initialize Input Bridge
  m_input_bridge = std::make_unique<OcctInputBridge>(
      this,
      myView,
      myContext,
      this,
      [this]() { this->queue_draw(); },
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

  // --- CONTEXT MENU CREATION (Right Click) ---

  // 1. Menu item action
  m_actionGroup = Gio::SimpleActionGroup::create();
  m_actionGroup->add_action("centre_graph", [this]() {
      if (m_contextMenuTargetId != -1) {
          signal_locate_in_tree.emit(m_contextMenuTargetId); // Send signal to the tree view
      }
  });
  insert_action_group("viewer", m_actionGroup);

  // 2. Menu appearance and items
  auto menuModel = Gio::Menu::create();
  menuModel->append("Locate in Tree (Centre Graph)", "viewer.centre_graph");
  m_contextMenu.set_menu_model(menuModel);
  m_contextMenu.set_parent(*this);
  m_contextMenu.set_has_arrow(false);

  // 3. Right-click detector
  m_rightClickGesture = Gtk::GestureClick::create();
  m_rightClickGesture->set_button(GDK_BUTTON_SECONDARY);
  m_rightClickGesture->signal_released().connect(sigc::mem_fun(*this, &OcctGtkGLAreaViewer::onRightClick));
  add_controller(m_rightClickGesture);
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

  std::cout << "--- KEYBOARD LOG ---" << std::endl;
  std::cout << "GTK KeyVal: " << theKeyVal << ", KeyCode: " << theKeyCode << std::endl;

  const Aspect_VKey aVKey = OcctGtkTools::gtkKey2VKey(theKeyVal, theKeyCode);

  std::cout << "OpenCASCADE VKey: " << aVKey << " (Escape is: " << Aspect_VKey_Escape << ")" << std::endl;

  if (aVKey == Aspect_VKey_Escape) {
      std::cout << ">>> ESCAPE DETECTED! Disabling clipping and clearing selection. <<<" << std::endl;
      disableClippingPlane();
  }

  if (aVKey == Aspect_VKey_UNKNOWN) {
      std::cout << "Unknown key, ignoring." << std::endl;
      return false;
  }

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
    const Gtk::TreeModelColumn<bool>& colVisible,
    std::map<int, TDF_Label>& labelMap,
    std::map<int, TopLoc_Location>& locMap,
    std::map<int, int>& depthMap,
    std::map<int, TopoDS_Shape>& subShapeMap,
    std::map<int, Handle(AIS_InteractiveObject)>& displayMap,
    std::map<int, std::vector<int>>& hierarchyMap,
    int parentId,
    int& idCounter,
    const TopLoc_Location& parentLoc,
    int currentTreeDepth)
{
    TopLoc_Location localLoc = shapeTool->GetLocation(label);
    TopLoc_Location globalLoc = parentLoc * localLoc;

    TDF_Label definition_label = label;
    if (shapeTool->IsReference(label)) shapeTool->GetReferredShape(label, definition_label);

    Handle(TDataStd_Name) aNodeName;
    Glib::ustring name_str = "Unnamed";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();
    else if (definition_label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();

    auto row = *store->append(children);
    row[colName] = name_str;

    TopoDS_Shape baseShape;
    shapeTool->GetShape(definition_label, baseShape);
    row[colType] = baseShape.IsNull() ? "Unknown" : shapeTypeToString(baseShape.ShapeType());

    int currentId = idCounter++;
    row[colId] = currentId;
    row[colVisible] = true;

    if (parentId != -1) hierarchyMap[parentId].push_back(currentId);
    labelMap[currentId] = definition_label;
    locMap[currentId]   = globalLoc;
    depthMap[currentId] = currentTreeDepth;

    TDF_LabelSequence comps;
    shapeTool->GetComponents(definition_label, comps);
    if (comps.Length() > 0) {
        for (Standard_Integer i = 1; i <= comps.Length(); i++) {
            processXCAFLabelRecursive(comps.Value(i), store, row.children(), shapeTool, colName, colType, colId, colVisible, labelMap, locMap, depthMap, subShapeMap, displayMap, hierarchyMap, currentId, idCounter, globalLoc, currentTreeDepth + 1);
        }
    }
    else {
        // --- 1. SHAPE ACQUISITION ---
        TopoDS_Shape partShape;
        shapeTool->GetShape(label, partShape);

        if (partShape.IsNull()) {
            return; // Geometry-less "ghost" node, skip
        }

        Handle(AIS_ColoredShape) leafPrs = new AIS_ColoredShape(partShape);
        leafPrs->SetLocalTransformation(parentLoc.Transformation());

        // --- 2. PRESENTATION SETTINGS ---
        Handle(Prs3d_Drawer) drawer = leafPrs->Attributes();
        drawer->SetFaceBoundaryDraw(true);
        drawer->SetFaceBoundaryAspect(new Prs3d_LineAspect(Quantity_NOC_BLACK, Aspect_TOL_SOLID, 1.0));
        Handle(Prs3d_ShadingAspect) aspect = drawer->ShadingAspect();
        if (aspect.IsNull()) { aspect = new Prs3d_ShadingAspect(); drawer->SetShadingAspect(aspect); }
        Graphic3d_MaterialAspect mat(Graphic3d_NOM_SHINY_PLASTIC);
        mat.SetSpecularColor(Quantity_NOC_WHITE);
        mat.SetShininess(0.9);
        aspect->SetMaterial(mat);

        // --- 3. COLOR HANDLING ---
        XCAFPrs_IndexedDataMapOfShapeStyle settings;
        XCAFPrs::CollectStyleSettings(label, TopLoc_Location(), settings);

        // Default body color
        Quantity_Color mainCol = Quantity_NOC_GRAY70;
        if (settings.Contains(partShape)) {
            const XCAFPrs_Style& mainStyle = settings.FindFromKey(partShape);
            if (mainStyle.IsSetColorSurf()) mainCol = mainStyle.GetColorSurf();
        }
        leafPrs->SetColor(mainCol);

        // Face-specific coloring
        for (Standard_Integer i = 1; i <= settings.Extent(); ++i) {
            const TopoDS_Shape& subShape = settings.FindKey(i);
            const XCAFPrs_Style& style = settings.FindFromIndex(i);

            if (style.IsSetColorSurf()) {
                leafPrs->SetCustomColor(subShape, style.GetColorSurf());
            }
        }

        displayMap[currentId] = leafPrs;

        // --- 4. TREE STORAGE FOR SOLIDS AND FACES ---
        if (!partShape.IsNull()) {
            TopTools_IndexedMapOfShape solidMap;
            TopExp::MapShapes(partShape, TopAbs_SOLID, solidMap);

            if (solidMap.Extent() > 1) {
                // MULTI-BODY PART
                int solidIdx = 1;
                for (int s = 1; s <= solidMap.Extent(); s++) {
                    const TopoDS_Shape& solid = solidMap(s);

                    auto solidRow = *store->append(row.children());
                    solidRow[colName] = Glib::ustring::compose("Solid %1", solidIdx++);
                    solidRow[colType] = "MULTIBODY";
                    int solidId = idCounter++;
                    solidRow[colId] = solidId;
                    solidRow[colVisible] = true;

                    subShapeMap[solidId] = solid;
                    hierarchyMap[currentId].push_back(solidId);

                    TopExp_Explorer exp(solid, TopAbs_FACE);
                    int faceIdx = 1;
                    for (; exp.More(); exp.Next()) {
                        const TopoDS_Shape& face = exp.Current();
                        auto faceRow = *store->append(solidRow.children());
                        faceRow[colName] = Glib::ustring::compose("Face %1", faceIdx++);
                        faceRow[colType] = "FACE";
                        int faceId = idCounter++;
                        faceRow[colId] = faceId;
                        faceRow[colVisible] = true;

                        subShapeMap[faceId] = face;
                        hierarchyMap[solidId].push_back(faceId);
                    }
                }
            } else {
                // SINGLE-BODY PART
                TopExp_Explorer exp(partShape, TopAbs_FACE);
                int faceIdx = 1;
                for (; exp.More(); exp.Next()) {
                    const TopoDS_Shape& face = exp.Current();
                    auto faceRow = *store->append(row.children());
                    faceRow[colName] = Glib::ustring::compose("Face %1", faceIdx++);
                    faceRow[colType] = "FACE";
                    int faceId = idCounter++;
                    faceRow[colId] = faceId;
                    faceRow[colVisible] = true;

                    subShapeMap[faceId] = face;
                    hierarchyMap[currentId].push_back(faceId);
                }
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

        myContext->Activate(obj, 0); // Whole part
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_FACE));   // Faces
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_EDGE));   // Edges
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_VERTEX)); // Vertices
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

    if (m_display_map.find(id) != m_display_map.end()) {
        myContext->AddOrRemoveSelected(m_display_map[id], false);
        foundInScene = true;
    }

    if (!foundInScene) {
        TopoDS_Shape shapeToSelect;
        gp_Trsf parentTransform;
        bool hasTransform = false;

        if (m_subshape_map.find(id) != m_subshape_map.end()) {
            shapeToSelect = m_subshape_map[id];

            int parentId = -1;
            for (auto const& [pId, children] : m_tree_hierarchy) {
                for (int childId : children) {
                    if (childId == id) { parentId = pId; break; }
                }
                if (parentId != -1) break;
            }

            if (parentId != -1 && m_display_map.find(parentId) != m_display_map.end()) {
                parentTransform = m_display_map[parentId]->LocalTransformation();
                hasTransform = true;
            }
        }

        if (!shapeToSelect.IsNull()) {
            m_selected_part_ais = new AIS_Shape(shapeToSelect);
            if (hasTransform) m_selected_part_ais->SetLocalTransformation(parentTransform);

            m_selected_part_ais->SetColor(Quantity_NOC_ORANGE);
            m_selected_part_ais->SetZLayer(Graphic3d_ZLayerId_Topmost);
            myContext->Display(m_selected_part_ais, false);
        }
    }
    myContext->UpdateCurrentViewer();
    queue_draw();
}

void OcctGtkGLAreaViewer::fitToShapeById(int id)
{
    TopoDS_Shape shapeToFit;
    if (m_label_map.count(id) && m_loc_map.count(id)) {
        Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());

        TDF_Label instanceLabel = m_label_map[id];
        TDF_Label defLabel = instanceLabel;
        if (shapeTool->IsReference(instanceLabel)) {
            shapeTool->GetReferredShape(instanceLabel, defLabel);
        }

        TopoDS_Shape rawShape;
        shapeTool->GetShape(defLabel, rawShape);
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

    // Direct memory match in subshapes
    for (auto const& [id, faceShape] : m_subshape_map) {
        if (selectedShape.IsSame(faceShape)) return id;
    }

    if (selectedShape.IsNull()) return -1;

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    int bestId = -1;
    int maxDepth = -1;
    for (auto const& [id, label] : m_label_map) {
        TopoDS_Shape partShape;
        shapeTool->GetShape(label, partShape);
        if (partShape.IsNull()) continue;

        if (partShape.IsSame(selectedShape)) {
            int currentDepth = 0;
            if (m_depth_map.find(id) != m_depth_map.end()) currentDepth = m_depth_map[id];
            if (currentDepth > maxDepth) { maxDepth = currentDepth; bestId = id; }
        }
    }
    return bestId;
}

// --- DIAGNOSTIC XCAF STRUCTURE DUMP ---
static void dumpRawXCAFStructure(const TDF_Label& label, const Handle(XCAFDoc_ShapeTool)& shapeTool, int depth = 0)
{
    std::string indent(depth * 4, ' ');

    TDF_Tool::Entry(label, entryStr);

    Handle(TDataStd_Name) aNodeName;
    std::string name = "Unnamed";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
        name = TCollection_AsciiString(aNodeName->Get()).ToCString();
    }

    std::string typeInfo = "";
    if (shapeTool->IsAssembly(label)) typeInfo += "[ASSEMBLY] ";
    if (shapeTool->IsSimpleShape(label)) typeInfo += "[SHAPE] ";

    if (shapeTool->IsReference(label)) {
        TDF_Label refLabel;
        shapeTool->GetReferredShape(label, refLabel);
        TCollection_AsciiString refEntryStr;
        TDF_Tool::Entry(refLabel, refEntryStr);
        typeInfo += "[REFERENCE (Instance) -> points to " + std::string(refEntryStr.ToCString()) + "] ";

        if (name == "Unnamed" && refLabel.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
            name = TCollection_AsciiString(aNodeName->Get()).ToCString() + std::string(" (Name from definition)");
        }
    }

    std::cout << indent << entryStr.ToCString() << " | " << name << " | " << typeInfo << std::endl;

    TDF_ChildIterator child_it(label, Standard_False);
    for (; child_it.More(); child_it.Next()) {
        dumpRawXCAFStructure(child_it.Value(), shapeTool, depth + 1);
    }
}

static void dumpColorStructure(const TDF_Label& label, const Handle(XCAFDoc_ShapeTool)& shapeTool, const Handle(XCAFDoc_ColorTool)& colorTool, int depth = 0)
{
    std::string indent(depth * 4, ' ');
    TCollection_AsciiString entryStr;
    TDF_Tool::Entry(label, entryStr);

    Handle(TDataStd_Name) aNodeName;
    std::string name = "Unnamed";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
        name = TCollection_AsciiString(aNodeName->Get()).ToCString();
    }

    TDF_Label defLabel = label;
    bool isInstance = shapeTool->IsReference(label);
    if (isInstance) {
        shapeTool->GetReferredShape(label, defLabel);
    }

    Quantity_Color col;
    std::string colorInfo = " | Color: NONE";

    if (colorTool->GetColor(label, XCAFDoc_ColorSurf, col) || colorTool->GetColor(label, XCAFDoc_ColorGen, col) ||
        colorTool->GetColor(defLabel, XCAFDoc_ColorSurf, col) || colorTool->GetColor(defLabel, XCAFDoc_ColorGen, col)) {
        char colorBuf[100];
        sprintf(colorBuf, " | Color: RGB(%.2f, %.2f, %.2f)", col.Red(), col.Green(), col.Blue());
        colorInfo = colorBuf;
    }

    std::string typeInfo = isInstance ? "[Instance -> searching definition]" : "[Assembly/Shape]";
    std::cout << indent << entryStr.ToCString() << " " << typeInfo << " " << name << colorInfo << std::endl;

    TDF_LabelSequence subLabels;
    shapeTool->GetSubShapes(defLabel, subLabels);
    if (subLabels.Length() > 0) {
        bool foundSubColor = false;
        for (Standard_Integer i = 1; i <= subLabels.Length(); i++) {
            TDF_Label subL = subLabels.Value(i);
            if (colorTool->GetColor(subL, XCAFDoc_ColorSurf, col) || colorTool->GetColor(subL, XCAFDoc_ColorGen, col)) {
                if (!foundSubColor) {
                    std::cout << indent << "  ---> Colored sub-labels (Faces) found:" << std::endl;
                    foundSubColor = true;
                }
                TCollection_AsciiString subEntry;
                TDF_Tool::Entry(subL, subEntry);
                std::cout << indent << "       " << subEntry.ToCString() << " | Color: RGB(" << col.Red() << ", " << col.Green() << ", " << col.Blue() << ")" << std::endl;
            }
        }
    }

    TDF_ChildIterator child_it(label, Standard_False);
    for (; child_it.More(); child_it.Next()) {
        dumpColorStructure(child_it.Value(), shapeTool, colorTool, depth + 1);
    }
}

bool OcctGtkGLAreaViewer::loadStepFile(const Glib::ustring& filePath,
                                Glib::RefPtr<Gtk::TreeStore>& treeModel,
                                const Gtk::TreeModelColumn<Glib::ustring>& colName,
                                const Gtk::TreeModelColumn<Glib::ustring>& colType,
                                const Gtk::TreeModelColumn<int>& colId,
                                const Gtk::TreeModelColumn<bool>& colVisible)
{
    auto updateStatus = [&](const std::string& msg) {
        signal_status_message.emit(msg);
        while (Glib::MainContext::get_default()->pending()) {
            Glib::MainContext::get_default()->iteration(false);
        }
    };

    updateStatus("Initializing document...");

    Handle(Prs3d_Drawer) globalDrawer = myContext->DefaultDrawer();
    globalDrawer->SetTypeOfDeflection(Aspect_TOD_RELATIVE);
    globalDrawer->SetDeviationCoefficient(0.001);
    globalDrawer->SetDeviationAngle(12.0 * M_PI / 180.0);

    globalDrawer->ShadingAspect()->SetColor(Quantity_NOC_GRAY70);
    Graphic3d_MaterialAspect globalMat(Graphic3d_NOM_SHINY_PLASTIC);
    globalMat.SetSpecularColor(Quantity_NOC_WHITE);
    globalMat.SetShininess(0.9);
    globalDrawer->ShadingAspect()->SetMaterial(globalMat);
    globalDrawer->SetFaceBoundaryDraw(true);
    globalDrawer->SetFaceBoundaryAspect(new Prs3d_LineAspect(Quantity_NOC_BLACK, Aspect_TOL_SOLID, 1.0));

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
        m_display_map.clear();
        m_tree_hierarchy.clear();
        m_next_id = 1;
        if(treeModel) treeModel->clear();
    }

    Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    TDF_LabelSequence rootsBefore;
    shapeTool->GetFreeShapes(rootsBefore);
    Standard_Integer countBefore = rootsBefore.Length();

    STEPCAFControl_Reader reader;
    reader.SetColorMode(true);
    reader.SetNameMode(true);
    reader.SetLayerMode(true);

    updateStatus("Phase 1/4: Reading STEP file from disk...");
    if (reader.ReadFile(path.ToCString()) != IFSelect_RetDone) {
        updateStatus("Error: Unable to read file!");
        return false;
    }

    updateStatus("Phase 2/4: Transferring geometry (0 %)...");
    Handle(GtkStepProgress) progress = new GtkStepProgress([&](int pct) {
        updateStatus("Phase 2/4: Transferring geometry (" + std::to_string(pct) + " %)");
    });

    if (!reader.Transfer(m_doc, progress->Start())) {
        updateStatus("Error: Geometry transfer failed!");
        return false;
    }

    std::cout << "\n================ START XCAF COLOR DIAGNOSTICS ================" << std::endl;
    Handle(XCAFDoc_ColorTool) colorTool = XCAFDoc_DocumentTool::ColorTool(m_doc->Main());
    TDF_LabelSequence freeShapes;
    shapeTool->GetFreeShapes(freeShapes);
    for (Standard_Integer i = 1; i <= freeShapes.Length(); i++) {
        dumpColorStructure(freeShapes.Value(i), shapeTool, colorTool, 0);
    }
    std::cout << "================ END XCAF COLOR DIAGNOSTICS ================\n" << std::endl;

    TDF_LabelSequence rootsAfter;
    shapeTool->GetFreeShapes(rootsAfter);

    updateStatus("Phase 3/4: Analyzing hierarchy...");
    std::string fileName = Glib::path_get_basename(filePath.c_str());

    auto fileRow = *treeModel->append();
    fileRow[colName] = fileName;
    fileRow[colType] = "FILE";
    int fileId = m_next_id++;
    fileRow[colId] = fileId;
    fileRow[colVisible] = true;

    for (Standard_Integer i = countBefore + 1; i <= rootsAfter.Length(); i++) {
        TDF_Label label = rootsAfter.Value(i);
        if(treeModel) {
            TopLoc_Location startLoc;
            processXCAFLabelRecursive(label, treeModel, fileRow.children(), shapeTool,
                                      colName, colType, colId, colVisible,
                                      m_label_map, m_loc_map, m_depth_map, m_subshape_map,
                                      m_display_map, m_tree_hierarchy, fileId,
                                      m_next_id, startLoc, 0);
        }
    }

    updateStatus("Phase 3.5/4: 3D Mesh Calculation (Tessellation)...");
    for (Standard_Integer i = countBefore + 1; i <= rootsAfter.Length(); i++) {
        TopoDS_Shape rootShape;
        shapeTool->GetShape(rootsAfter.Value(i), rootShape);
        if (!rootShape.IsNull()) {
            BRepMesh_IncrementalMesh mesher(rootShape, 0.001, Standard_True, 12.0 * M_PI / 180.0, Standard_True);
        }
    }

    updateStatus("Phase 4/4: Sending to GPU and Creating Selection...");
    int totalNewItems = 0;
    for (auto const& [id, aisObj] : m_display_map) {
        if (id >= fileId) totalNewItems++;
    }

    int currentItem = 0;
    for (auto const& [id, aisObj] : m_display_map) {
        if (id < fileId) continue;
        currentItem++;

        if (currentItem % 10 == 0 || currentItem == totalNewItems || currentItem == 1) {
            int percent = (currentItem * 100) / (totalNewItems > 0 ? totalNewItems : 1);
            updateStatus("Phase 4/4: GPU & Selection (" + std::to_string(percent) + " %) - Part " + std::to_string(currentItem) + " of " + std::to_string(totalNewItems));
        }

        myContext->Display(aisObj, AIS_Shaded, 0, false);
        myContext->Activate(aisObj, 0);
        myContext->Activate(aisObj, AIS_Shape::SelectionMode(TopAbs_FACE));
    }

    myContext->UpdateCurrentViewer();
    myView->FitAll();
    drawOriginAxes();
    updateStatus("Model loaded successfully.");
    return true;
}

void OcctGtkGLAreaViewer::toggleClippingPlane()
{
    if (myView.IsNull() || myContext.IsNull()) return;

    // --- DOUBLE CLICK DETECTION ---
    static gint64 last_click_time = 0;
    gint64 current_time = g_get_monotonic_time();
    bool isDoubleClick = (current_time - last_click_time < 400000); // 400ms threshold
    last_click_time = isDoubleClick ? 0 : current_time;

    // 1. Initial plane creation
    if (m_clipPlane.IsNull()) {
        m_clipPlane = new Graphic3d_ClipPlane(gp_Pln(gp_Pnt(0,0,0), gp_Dir(1,0,0)));
        m_clipPlane->SetCapping(true);

        Quantity_Color darkRed(0.55, 0.0, 0.0, Quantity_TOC_RGB);
        Handle(Graphic3d_AspectFillArea3d) cappingAspect = m_clipPlane->CappingAspect();
        cappingAspect->SetInteriorColor(darkRed);

        Graphic3d_MaterialAspect mat = cappingAspect->FrontMaterial();
        mat.SetColor(darkRed);
        mat.SetAmbientColor(darkRed);
        cappingAspect->SetFrontMaterial(mat);

        myView->AddClipPlane(m_clipPlane);
    }

    // --- DOUBLE CLICK ACTION (Flip Plane) ---
    if (isDoubleClick) {
        gp_Pln currentPln = m_clipPlane->ToPlane();
        gp_Ax1 axis = currentPln.Axis();
        axis.SetDirection(axis.Direction().Reversed());
        currentPln.SetAxis(axis);
        m_clipPlane->SetEquation(currentPln);
        m_clipPlane->SetOn(true);
        myView->Redraw();
        return;
    }

    // --- SINGLE CLICK ACTION ---
    gp_Dir finalDir(0, 0, 1);
    gp_Pnt finalPos(0, 0, 0);
    bool hasExplicitDir = false;
    bool hasExplicitPos = false;
    gp_Dir fallbackDir(0, 0, 1);
    bool hasFallbackDir = false;

    std::vector<gp_Pnt> selectedVertices;

    for (myContext->InitSelected(); myContext->MoreSelected(); myContext->NextSelected())
    {
        if (!myContext->HasSelectedShape()) continue;
        TopoDS_Shape shape = myContext->SelectedShape();

        if (shape.ShapeType() == TopAbs_FACE) {
            TopoDS_Face face = TopoDS::Face(shape);
            BRepAdaptor_Surface surf(face);
            if (surf.GetType() == GeomAbs_Plane) {
                gp_Pln pln = surf.Plane();
                finalDir = pln.Axis().Direction();
                hasExplicitDir = true;
                if (!hasExplicitPos) finalPos = pln.Location();
            } else if (surf.GetType() == GeomAbs_Cylinder) {
                gp_Cylinder cyl = surf.Cylinder();
                finalPos = cyl.Location();
                hasExplicitPos = true;
                fallbackDir = cyl.Axis().Direction();
                hasFallbackDir = true;
            }
        }
        else if (shape.ShapeType() == TopAbs_EDGE) {
            TopoDS_Edge edge = TopoDS::Edge(shape);
            BRepAdaptor_Curve curve(edge);
            if (curve.GetType() == GeomAbs_Line) {
                gp_Lin lin = curve.Line();
                finalDir = lin.Direction();
                hasExplicitDir = true;
                if (!hasExplicitPos) finalPos = lin.Location();
            } else if (curve.GetType() == GeomAbs_Circle) {
                gp_Circ circ = curve.Circle();
                finalPos = circ.Location();
                hasExplicitPos = true;
                fallbackDir = circ.Axis().Direction();
                hasFallbackDir = true;
            }
        }
        else if (shape.ShapeType() == TopAbs_VERTEX) {
            TopoDS_Vertex vertex = TopoDS::Vertex(shape);
            gp_Pnt pnt = BRep_Tool::Pnt(vertex);
            selectedVertices.push_back(pnt);
            finalPos = pnt;
            hasExplicitPos = true;
        }
    }

    // --- 3-POINT PLANE LOGIC ---
    if (selectedVertices.size() == 3) {
        gp_Vec v1(selectedVertices[0], selectedVertices[1]);
        gp_Vec v2(selectedVertices[0], selectedVertices[2]);
        gp_Vec normal = v1.Crossed(v2);

        if (normal.Magnitude() > gp::Resolution()) {
            finalDir = gp_Dir(normal);
            finalPos = selectedVertices[0];
            hasExplicitDir = true;
            hasExplicitPos = true;
            std::cout << "[CLIP] Plane created from 3 points." << std::endl;
        } else {
            std::cout << "[CLIP] ERROR: Selected 3 points are collinear!" << std::endl;
        }
    }
    else if (!hasExplicitDir && hasFallbackDir) {
        finalDir = fallbackDir;
        hasExplicitDir = true;
    }

    if (hasExplicitDir || hasExplicitPos) {
        Standard_Real vx, vy, vz;
        myView->Proj(vx, vy, vz);
        gp_Dir viewDir(vx, vy, vz);

        if (finalDir.Dot(viewDir) > 0.0) {
            finalDir.Reverse();
        }

        gp_Pln newPlane(finalPos, finalDir);
        m_clipPlane->SetEquation(newPlane);
        m_clipPlane->SetOn(true);
        myContext->ClearSelected(false);

    } else {
        // Toggle behavior if nothing selected
        bool isCurrentlyOn = m_clipPlane->IsOn();
        m_clipPlane->SetOn(!isCurrentlyOn);
    }

    myView->Redraw();
}

void OcctGtkGLAreaViewer::disableClippingPlane()
{
    if (!m_clipPlane.IsNull()) {
        m_clipPlane->SetOn(false);
        m_clipState = 0;
    }
    if (!myContext.IsNull()) {
        myContext->ClearSelected(false);
    }
    if (!myView.IsNull()) {
        myView->Redraw();
    }
}

void OcctGtkGLAreaViewer::setShapeVisibility(int id, bool visible)
{
    if (myContext.IsNull()) return;

    bool screenNeedsUpdate = false;
    bool isMajorHide = false;

    // 1. Hiding whole main parts
    if (m_display_map.find(id) != m_display_map.end()) {
        Handle(AIS_InteractiveObject) obj = m_display_map[id];
        if (!obj.IsNull()) {
            if (visible) {
                myContext->Display(obj, AIS_Shaded, 0, false);
                myContext->Activate(obj, 0);
                myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_EDGE));
                myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_VERTEX));
                myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_FACE));
            } else {
                myContext->Erase(obj, false);
            }
            screenNeedsUpdate = true;
            isMajorHide = true;
        }
    }
    // 2. Hiding sub-shapes (Multibody solids or faces)
    else if (m_subshape_map.find(id) != m_subshape_map.end()) {
        TopoDS_Shape subShape = m_subshape_map[id];

        int rootDisplayId = -1;
        int currentSearchId = id;

        while (rootDisplayId == -1) {
            int currentParent = -1;
            for (auto const& [pId, children] : m_tree_hierarchy) {
                for (int childId : children) {
                    if (childId == currentSearchId) { currentParent = pId; break; }
                }
                if (currentParent != -1) break;
            }
            if (currentParent == -1) break;
            if (m_display_map.find(currentParent) != m_display_map.end()) {
                rootDisplayId = currentParent;
            } else {
                currentSearchId = currentParent;
            }
        }

        if (rootDisplayId != -1) {
            Handle(AIS_ColoredShape) coloredObj = Handle(AIS_ColoredShape)::DownCast(m_display_map[rootDisplayId]);
            if (!coloredObj.IsNull() && !coloredObj->Shape().IsNull()) {

                TopExp_Explorer expFacesToHide(subShape, TopAbs_FACE);
                bool madeChanges = false;

                for (; expFacesToHide.More(); expFacesToHide.Next()) {
                    const TopoDS_Shape& faceToHide = expFacesToHide.Current();
                    TopExp_Explorer expAis(coloredObj->Shape(), TopAbs_FACE);
                    for (; expAis.More(); expAis.Next()) {
                        if (expAis.Current().IsPartner(faceToHide)) {
                            coloredObj->CustomAspects(expAis.Current())->SetHidden(!visible);
                            madeChanges = true;
                            break;
                        }
                    }
                }

                if (madeChanges) {
                    myContext->Redisplay(coloredObj, Standard_True);
                    screenNeedsUpdate = true;
                }

                if (subShape.ShapeType() == TopAbs_SOLID) {
                    isMajorHide = true;
                }
            }
        }
    }

    // 3. TREE RECURSION
    if (m_tree_hierarchy.find(id) != m_tree_hierarchy.end()) {
        for (int childId : m_tree_hierarchy[id]) {
            if (isMajorHide && m_subshape_map.find(childId) != m_subshape_map.end()) {
                continue;
            }
            setShapeVisibility(childId, visible);
        }
    }

    if (screenNeedsUpdate) {
        myContext->UpdateCurrentViewer();
        queue_draw();
    }
}

void OcctGtkGLAreaViewer::drawOriginAxes()
{
    if (myContext.IsNull()) return;

    gp_Pnt origin(0, 0, 0);
    Handle(AIS_Shape) aisX = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(30, 0, 0)));
    Handle(AIS_Shape) aisY = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(0, 30, 0)));
    Handle(AIS_Shape) aisZ = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(0, 0, 30)));

    aisX->SetMutable(false);
    aisY->SetMutable(false);
    aisZ->SetMutable(false);

    aisX->SetColor(Quantity_NOC_RED);
    aisY->SetColor(Quantity_NOC_GREEN);
    aisZ->SetColor(Quantity_NOC_BLUE);

    aisX->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_RED, Aspect_TOL_SOLID, 2.0));
    aisY->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_GREEN, Aspect_TOL_SOLID, 2.0));
    aisZ->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_BLUE, Aspect_TOL_SOLID, 2.0));

    aisX->SetZLayer(Graphic3d_ZLayerId_Topmost);
    aisY->SetZLayer(Graphic3d_ZLayerId_Topmost);
    aisZ->SetZLayer(Graphic3d_ZLayerId_Topmost);

    myContext->Display(aisX, 0, -1, false);
    myContext->Display(aisY, 0, -1, false);
    myContext->Display(aisZ, 0, -1, false);
}

void OcctGtkGLAreaViewer::onRightClick(int n_press, double x, double y)
{
    int selectedId = getSelectedShapeId();

    if (selectedId != -1) {
        m_contextMenuTargetId = selectedId;

        Gdk::Rectangle rect;
        rect.set_x(static_cast<int>(x));
        rect.set_y(static_cast<int>(y));
        rect.set_width(1);
        rect.set_height(1);

        m_contextMenu.set_pointing_to(rect);
        m_contextMenu.popup();
    }
}
