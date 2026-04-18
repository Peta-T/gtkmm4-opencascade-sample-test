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
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <glibmm/main.h>
#include <Geom_Axis2Placement.hxx>
#include <AIS_Trihedron.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TDF_Tool.hxx>
#include <iostream>
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


/// #include <XCAFPrs_DataMapOfShapeStyle.hxx>

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

// Pomocná třída pro sledování průběhu OpenCASCADE operací
class GtkStepProgress : public Message_ProgressIndicator {
private:
    std::function<void(int)> m_callback;
    int m_lastPercent = -1;
public:
    GtkStepProgress(std::function<void(int)> cb) : m_callback(cb) {}

protected:
    virtual void Show(const Message_ProgressScope& theScope, const Standard_Boolean isForce) override {
        // GetPosition() vrací hodnotu 0.0 až 1.0
        int percent = static_cast<int>(GetPosition() * 100.0);

        // Aktualizujeme UI pouze když se procento změní (abychom nebrzdili výpočet)
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

  Handle(Prs3d_Drawer) ctxDrawer = myContext->DefaultDrawer();

  // Odchylka bude relativní k velikosti daného dílu
  ctxDrawer->SetTypeOfDeflection(Aspect_TOD_RELATIVE);

  // Vzdálenostní odchylka (výchozí bývá 0.001). Čím menší, tím přesnější a detailnější.
  ctxDrawer->SetDeviationCoefficient(0.001);

  // Úhlová odchylka pro zakřivené plochy (např. 10 stupňů místo výchozích 20).
  ctxDrawer->SetDeviationAngle(12.0*M_PI / 180.0);
  // Parametry: tvar, lineární deflekce, isRelative, úhlová deflekce, paralelní zpracování
///BRepMesh_IncrementalMesh mesher(rootShape, 0.001, Standard_True, 12.0 * M_PI / 180.0, Standard_True);

myContext->HighlightStyle(Prs3d_TypeOfHighlight_Dynamic)->SetColor(Quantity_NOC_ORANGE);
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_LocalDynamic)->SetColor(Quantity_NOC_ORANGE);

  // 2. Barva po kliknutí/výběru (můžeme dát tmavší oranžovou nebo nechat stejnou)
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_Selected)->SetColor(Quantity_NOC_ORANGE);
  myContext->HighlightStyle(Prs3d_TypeOfHighlight_LocalSelected)->SetColor(Quantity_NOC_ORANGE);

  myViewCube = new AIS_ViewCube();
  myViewCube->SetSize(50.0);        // Výchozí velikost bývá kolem 55-60. Zkuste 30 až 40.
    myViewCube->SetFontHeight(10.0);  // Zmenšíme text, aby se vešel na menší stěny (výchozí je 16).
   // myViewCube->SetAxesPadding(2.0);  // Zmenšíme mezeru mezi krychlí a šipkami.

    // Volitelně: Pokud byste chtěl změnit i průhlednost (0.0 = neprůhledná, 1.0 = zcela průhledná)
    myViewCube->SetTransparency(0.8);

  myViewCube->SetViewAnimation(myViewAnimation);
  myViewCube->SetFixedAnimationLoop(false);
  myViewCube->SetAutoStartAnimation(true);

   Handle(Prs3d_DatumAspect) datumAspect = myViewCube->Attributes()->DatumAspect();
  if (datumAspect.IsNull()) {
        datumAspect = new Prs3d_DatumAspect();
        myViewCube->Attributes()->SetDatumAspect(datumAspect);
    }

    // Nastavení černé barvy pro text všech tří os
    datumAspect->TextAspect(Prs3d_DatumParts_XAxis)->SetColor(Quantity_NOC_RED);
    datumAspect->TextAspect(Prs3d_DatumParts_YAxis)->SetColor(Quantity_NOC_GREEN);
    datumAspect->TextAspect(Prs3d_DatumParts_ZAxis)->SetColor(Quantity_NOC_BLUE);


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

  // --- VYTVOŘENÍ KONTEXTOVÉHO MENU (Pravé tlačítko) ---

  // 1. Akce, která se stane po kliknutí na položku v menu
  m_actionGroup = Gio::SimpleActionGroup::create();
  m_actionGroup->add_action("centre_graph", [this]() {
      if (m_contextMenuTargetId != -1) {
          signal_locate_in_tree.emit(m_contextMenuTargetId); // Pošleme signál dál!
      }
  });
  insert_action_group("viewer", m_actionGroup);

  // 2. Vzhled a položky menu
  auto menuModel = Gio::Menu::create();
  menuModel->append("Najít ve stromu (Centre Graph)", "viewer.centre_graph");
  m_contextMenu.set_menu_model(menuModel);
  m_contextMenu.set_parent(*this);
  m_contextMenu.set_has_arrow(false); // Vypadá to víc jako klasické Win/Mac menu

  // 3. Detektor pravého kliknutí (GDK_BUTTON_SECONDARY)
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

  // 1. LOG: Jaký surový kód nám poslal GTK?
  std::cout << "--- LOG KLÁVESNICE ---" << std::endl;
  std::cout << "GTK KeyVal: " << theKeyVal << ", KeyCode: " << theKeyCode << std::endl;

  const Aspect_VKey aVKey = OcctGtkTools::gtkKey2VKey(theKeyVal, theKeyCode);

  // 2. LOG: Jak si to přeložil OpenCASCADE?
  std::cout << "OpenCASCADE VKey: " << aVKey << " (Escape má číslo: " << Aspect_VKey_Escape << ")" << std::endl;

  if (aVKey == Aspect_VKey_Escape) {
      // 3. LOG: Úspěšně jsme chytili ESC
      std::cout << ">>> DETEKOVÁN ESCAPE! Vypínám řez a mažu výběr. <<<" << std::endl;
      disableClippingPlane();
  }

  if (aVKey == Aspect_VKey_UNKNOWN) {
      std::cout << "Neznámá klávesa, ignoruji." << std::endl;
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
    TopLoc_Location myLocalLoc = shapeTool->GetLocation(label);
    TopLoc_Location myGlobalLoc = parentLoc * myLocalLoc;

    TDF_Label definition_label = label;
    if (shapeTool->IsReference(label)) shapeTool->GetReferredShape(label, definition_label);

    Handle(TDataStd_Name) aNodeName;
    Glib::ustring name_str = "Unnamed";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();
    else if (definition_label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) name_str = TCollection_AsciiString(aNodeName->Get()).ToCString();

    auto row = *store->append(children);
    row[colName] = name_str;

    TopoDS_Shape baseShape; // Čistý tvar bez posunu
    shapeTool->GetShape(definition_label, baseShape);
    row[colType] = baseShape.IsNull() ? "Unknown" : shapeTypeToString(baseShape.ShapeType());

    int currentId = idCounter++;
    row[colId] = currentId;
    row[colVisible] = true;

    if (parentId != -1) hierarchyMap[parentId].push_back(currentId);
    labelMap[currentId] = definition_label;
    locMap[currentId]   = myGlobalLoc;
    depthMap[currentId] = currentTreeDepth;

    TDF_LabelSequence comps;
    shapeTool->GetComponents(definition_label, comps);
    if (comps.Length() > 0) {
        for (Standard_Integer i = 1; i <= comps.Length(); i++) {
            processXCAFLabelRecursive(comps.Value(i), store, row.children(), shapeTool, colName, colType, colId, colVisible, labelMap, locMap, depthMap, subShapeMap, displayMap, hierarchyMap, currentId, idCounter, myGlobalLoc, currentTreeDepth + 1);
        }
    }
    else {
        // --- 1. ZÍSKÁNÍ TVARU ---
        TopoDS_Shape partShape;
        shapeTool->GetShape(label, partShape);

        if (partShape.IsNull()) {
            return; // Je to jen "duch" bez geometrie, zahodíme ho a jdeme dál
        }

        // Náš poslušný objekt, který umí skrývat plošky
        Handle(AIS_ColoredShape) leafPrs = new AIS_ColoredShape(partShape);
        leafPrs->SetLocalTransformation(parentLoc.Transformation());

        // --- 2. KOSMETIKA ---
        Handle(Prs3d_Drawer) drawer = leafPrs->Attributes();
        drawer->SetFaceBoundaryDraw(true);
        drawer->SetFaceBoundaryAspect(new Prs3d_LineAspect(Quantity_NOC_BLACK, Aspect_TOL_SOLID, 1.0));
        Handle(Prs3d_ShadingAspect) aspect = drawer->ShadingAspect();
        if (aspect.IsNull()) { aspect = new Prs3d_ShadingAspect(); drawer->SetShadingAspect(aspect); }
        Graphic3d_MaterialAspect mat(Graphic3d_NOM_SHINY_PLASTIC);
        mat.SetSpecularColor(Quantity_NOC_WHITE);
        mat.SetShininess(0.9);
        aspect->SetMaterial(mat);

        // --- 3. BAREVNÁ MAGIE - PŘESNÝ KLON NATIVNÍHO CHOVÁNÍ OCCT ---
        XCAFPrs_IndexedDataMapOfShapeStyle settings;
        XCAFPrs::CollectStyleSettings(label, TopLoc_Location(), settings);

        // Výchozí barva tělesa
        Quantity_Color mainCol = Quantity_NOC_GRAY70;
        if (settings.Contains(partShape)) {
            const XCAFPrs_Style& mainStyle = settings.FindFromKey(partShape);
            if (mainStyle.IsSetColorSurf()) mainCol = mainStyle.GetColorSurf();
        }
        leafPrs->SetColor(mainCol);

        // Obarvení plošek! Zde používáme nativní metodu SetCustomColor.
        // Tato metoda přesně propojila FreeCAD a XCAFPrs_AISObject.
        for (Standard_Integer i = 1; i <= settings.Extent(); ++i) {
            const TopoDS_Shape& subShape = settings.FindKey(i);
            const XCAFPrs_Style& style = settings.FindFromIndex(i);

            if (style.IsSetColorSurf()) {
                leafPrs->SetCustomColor(subShape, style.GetColorSurf());
            }
        }

        displayMap[currentId] = leafPrs;

        // --- 4. CHYTRÉ ULOŽENÍ TĚLES A PLOCH PRO STROM ---
        if (!partShape.IsNull()) {
            TopTools_IndexedMapOfShape solidMap;
            TopExp::MapShapes(partShape, TopAbs_SOLID, solidMap);

            if (solidMap.Extent() > 1) {
                // MULTI-BODY PART
                int solidIdx = 1;
                for (int s = 1; s <= solidMap.Extent(); s++) {
                    const TopoDS_Shape& solid = solidMap(s);

                    // Čisté pojmenování a typ MULTIBODY
                    auto solidRow = *store->append(row.children());
                    solidRow[colName] = Glib::ustring::compose("Těleso %1", solidIdx++);
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
                        faceRow[colName] = Glib::ustring::compose("Plocha %1", faceIdx++);
                        faceRow[colType] = "FACE";
                        int faceId = idCounter++;
                        faceRow[colId] = faceId;
                        faceRow[colVisible] = true;

                        subShapeMap[faceId] = face;
                        hierarchyMap[solidId].push_back(faceId);
                    }
                }
            } else {
                // KLASICKÝ JEDNODÍL
                TopExp_Explorer exp(partShape, TopAbs_FACE);
                int faceIdx = 1;
                for (; exp.More(); exp.Next()) {
                    const TopoDS_Shape& face = exp.Current();
                    auto faceRow = *store->append(row.children());
                    faceRow[colName] = Glib::ustring::compose("Plocha %1", faceIdx++);
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

        // Záměrně už neděláme Deactivate(obj)!
        // Necháme aktivní výběr celého dílu (mód 0) a přidáme k němu sub-módy.
        myContext->Activate(obj, 0); // Celý díl
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_FACE));   // Plochy
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_EDGE));   // Hrany
        myContext->Activate(obj, AIS_Shape::SelectionMode(TopAbs_VERTEX)); // Body
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
            // K lokaci z instance ještě přidáme virtuální posun rodiče
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
        Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());

        // --- OPRAVA DVOJITÉ TRANSFORMACE PRO ZOOM ---
        TDF_Label instanceLabel = m_label_map[id];
        TDF_Label defLabel = instanceLabel;
        if (ShapeTool->IsReference(instanceLabel)) {
            ShapeTool->GetReferredShape(instanceLabel, defLabel);
        }

        TopoDS_Shape rawShape;
        ShapeTool->GetShape(defLabel, rawShape);
        if (!rawShape.IsNull()) {
            shapeToFit = rawShape.Moved(m_loc_map[id]);
        }
        // ---------------------------------------------
    }
    else if (m_subshape_map.count(id)) {
        shapeToFit = m_subshape_map[id];
    }
    if (shapeToFit.IsNull()) return;

    // Tady musíme zahrnout #include <Bnd_Box.hxx> a #include <BRepBndLib.hxx>
    // (Předpokládám, že je tam nahoře už máte)
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

    // Přímá paměťová shoda v plochách (Nyní bude fungovat na 100 %)
    for (auto const& [id, faceShape] : m_subshape_map) {
        if (selectedShape.IsSame(faceShape)) return id;
    }

    if (selectedShape.IsNull()) return -1;

    Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    int bestId = -1;
    int maxDepth = -1;
    for (auto const& [id, label] : m_label_map) {
        TopoDS_Shape partShape;
        ShapeTool->GetShape(label, partShape);
        if (partShape.IsNull()) continue;

        if (partShape.IsSame(selectedShape)) {
            int currentDepth = 0;
            if (m_depth_map.find(id) != m_depth_map.end()) currentDepth = m_depth_map[id];
            if (currentDepth > maxDepth) { maxDepth = currentDepth; bestId = id; }
        }
    }
    return bestId;
}

// --- DIAGNOSTICKÝ RENTGEN XCAF STRUKTURY ---
static void dumpRawXCAFStructure(const TDF_Label& label, const Handle(XCAFDoc_ShapeTool)& shapeTool, int depth = 0)
{
    std::string indent(depth * 4, ' '); // Odsazení pro vizuální strom

    // 1. Získáme vnitřní adresu (např. 0:1:1:1)
    TCollection_AsciiString entryStr;
    TDF_Tool::Entry(label, entryStr);

    // 2. Získáme název
    Handle(TDataStd_Name) aNodeName;
    std::string name = "Bez_Jmena";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
        name = TCollection_AsciiString(aNodeName->Get()).ToCString();
    }

    // 3. Zjistíme typ z pohledu OpenCASCADE
    std::string typeInfo = "";
    if (shapeTool->IsAssembly(label)) typeInfo += "[SESTAVA] ";
    if (shapeTool->IsSimpleShape(label)) typeInfo += "[TVAR] ";

    if (shapeTool->IsReference(label)) {
        TDF_Label refLabel;
        shapeTool->GetReferredShape(label, refLabel);
        TCollection_AsciiString refEntryStr;
        TDF_Tool::Entry(refLabel, refEntryStr);
        typeInfo += "[ODKAZ (Instance) -> ukazuje na " + std::string(refEntryStr.ToCString()) + "] ";

        // Pokud je to odkaz a nemá vlastní jméno, zkusíme vzít jméno originálu
        if (name == "Bez_Jmena" && refLabel.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
            name = TCollection_AsciiString(aNodeName->Get()).ToCString() + std::string(" (Jméno z originálu)");
        }
    }

    // Vypíšeme řádek do terminálu
    std::cout << indent << entryStr.ToCString() << " | " << name << " | " << typeInfo << std::endl;

    // 4. Rekurzivně projdeme děti tohoto uzlu
    TDF_ChildIterator child_it(label, Standard_False);
    for (; child_it.More(); child_it.Next()) {
        dumpRawXCAFStructure(child_it.Value(), shapeTool, depth + 1);
    }
}
// -------------------------------------------

static void dumpColorStructure(const TDF_Label& label, const Handle(XCAFDoc_ShapeTool)& shapeTool, const Handle(XCAFDoc_ColorTool)& colorTool, int depth = 0)
{
    std::string indent(depth * 4, ' ');

    TCollection_AsciiString entryStr;
    TDF_Tool::Entry(label, entryStr);

    Handle(TDataStd_Name) aNodeName;
    std::string name = "Bez_Jmena";
    if (label.FindAttribute(TDataStd_Name::GetID(), aNodeName)) {
        name = TCollection_AsciiString(aNodeName->Get()).ToCString();
    }

    // --- KLÍČOVÁ ZMĚNA: Zjistíme, zda je to odkaz, a najdeme Originál ---
    TDF_Label defLabel = label;
    bool isInstance = shapeTool->IsReference(label);
    if (isInstance) {
        shapeTool->GetReferredShape(label, defLabel);
    }

    Quantity_Color col;
    std::string colorInfo = " | Barva: ZADNA";

    // Zkusíme najít barvu na instanci, a pokud tam není, hledáme na Originálu!
    if (colorTool->GetColor(label, XCAFDoc_ColorSurf, col) || colorTool->GetColor(label, XCAFDoc_ColorGen, col) ||
        colorTool->GetColor(defLabel, XCAFDoc_ColorSurf, col) || colorTool->GetColor(defLabel, XCAFDoc_ColorGen, col)) {
        char colorBuf[100];
        sprintf(colorBuf, " | Barva: RGB(%.2f, %.2f, %.2f)", col.Red(), col.Green(), col.Blue());
        colorInfo = colorBuf;
    }

    std::string typeInfo = isInstance ? "[Instance -> prohledávám originál]" : "[Sestava/Tvar]";
    std::cout << indent << entryStr.ToCString() << " " << typeInfo << " " << name << colorInfo << std::endl;

    // --- KLÍČOVÁ ZMĚNA 2: Prohledáme i pod-štítky originálu (jednotlivé plošky) ---
    TDF_LabelSequence subLabels;
    shapeTool->GetSubShapes(defLabel, subLabels);
    if (subLabels.Length() > 0) {
        bool foundSubColor = false;
        for (Standard_Integer i = 1; i <= subLabels.Length(); i++) {
            TDF_Label subL = subLabels.Value(i);
            if (colorTool->GetColor(subL, XCAFDoc_ColorSurf, col) || colorTool->GetColor(subL, XCAFDoc_ColorGen, col)) {
                if (!foundSubColor) {
                    std::cout << indent << "  ---> Nalezeny obarvené pod-štítky (Plošky):" << std::endl;
                    foundSubColor = true;
                }
                TCollection_AsciiString subEntry;
                TDF_Tool::Entry(subL, subEntry);
                std::cout << indent << "       " << subEntry.ToCString() << " | Barva: RGB(" << col.Red() << ", " << col.Green() << ", " << col.Blue() << ")" << std::endl;
            }
        }
    }

    // Rekurze do dětí (pokud je to sestava)
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

    updateStatus("Inicializace dokumentu...");

    // --- GLOBÁLNÍ VZHLED A KVALITA SÍTĚ ---
    Handle(Prs3d_Drawer) globalDrawer = myContext->DefaultDrawer();

    // Zásadní pro výkon: Relativní odchylka a úhel 12 stupňů
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

    Handle(XCAFDoc_ShapeTool) ShapeTool = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    TDF_LabelSequence rootsBefore;
    ShapeTool->GetFreeShapes(rootsBefore);
    Standard_Integer countBefore = rootsBefore.Length();

    STEPCAFControl_Reader reader;
    reader.SetColorMode(true);
    reader.SetNameMode(true);
    reader.SetLayerMode(true);

    updateStatus("Fáze 1/4: Čtení STEP souboru z disku...");
    if (reader.ReadFile(path.ToCString()) != IFSelect_RetDone) {
        updateStatus("Chyba: Soubor nelze přečíst!");
        return false;
    }

    updateStatus("Fáze 2/4: Převod geometrie (0 %)...");

    // Vytvoříme náš odposlouchávač a řekneme mu, ať volá updateStatus
    Handle(GtkStepProgress) progress = new GtkStepProgress([&](int pct) {
        updateStatus("Fáze 2/4: Převod geometrie (" + std::to_string(pct) + " %)");
    });

    // Spustíme Transfer a předáme mu náš progress (přes metodu Start())
    if (!reader.Transfer(m_doc, progress->Start())) {
        updateStatus("Chyba: Převod geometrie selhal!");
        return false;
    }

    // --- BAREVNÁ DIAGNOSTIKA (Volitelné, můžete smazat, pokud log nepotřebujete) ---
    std::cout << "\n================ START XCAF BAREVNE DIAGNOSTIKY ================" << std::endl;
    Handle(XCAFDoc_ColorTool) ColorTool = XCAFDoc_DocumentTool::ColorTool(m_doc->Main());
    TDF_LabelSequence freeShapes;
    ShapeTool->GetFreeShapes(freeShapes);
    for (Standard_Integer i = 1; i <= freeShapes.Length(); i++) {
        dumpColorStructure(freeShapes.Value(i), ShapeTool, ColorTool, 0);
    }
    std::cout << "================ KONEC XCAF BAREVNE DIAGNOSTIKY ================\n" << std::endl;
    // ---------------------------------------------------------------------------------

    TDF_LabelSequence rootsAfter;
    ShapeTool->GetFreeShapes(rootsAfter);

    updateStatus("Fáze 3/4: Analýza hierarchie...");

    // 1. Získání čistého názvu souboru (bezpečně přes Glib)
    std::string fileName = Glib::path_get_basename(filePath.c_str());

    // 2. Vytvoření kořenové položky FILE
    auto fileRow = *treeModel->append();
    fileRow[colName] = fileName;
    fileRow[colType] = "FILE";
    int fileId = m_next_id++;
    fileRow[colId] = fileId;
    fileRow[colVisible] = true;

    // 3. Rekurzivní analýza - směřujeme pod fileRow
    for (Standard_Integer i = countBefore + 1; i <= rootsAfter.Length(); i++) {
        TDF_Label label = rootsAfter.Value(i);
        if(treeModel) {
            TopLoc_Location startLoc;
            processXCAFLabelRecursive(label, treeModel, fileRow.children(), ShapeTool,
                                      colName, colType, colId, colVisible,
                                      m_label_map, m_loc_map, m_depth_map, m_subshape_map,
                                      m_display_map, m_tree_hierarchy, fileId,
                                      m_next_id, startLoc, 0);
        }
    }

    // --- FÁZE 3.5: VÝPOČET 3D SÍTĚ S RELATIVNÍ ODCHYLKOU ---
    updateStatus("Fáze 3.5/4: Výpočet 3D sítě (Tessellation)...");
    for (Standard_Integer i = countBefore + 1; i <= rootsAfter.Length(); i++) {
        TopoDS_Shape rootShape;
        ShapeTool->GetShape(rootsAfter.Value(i), rootShape);
        if (!rootShape.IsNull()) {
            // Relativní výpočet: Standard_True u druhého boolean parametru!
            BRepMesh_IncrementalMesh mesher(rootShape, 0.001, Standard_True, 12.0 * M_PI / 180.0, Standard_True);
        }
    }

    // --- FÁZE 4: ODESÍLÁNÍ NA GPU A TVORBA VÝBĚRU (V JEDNOM KROKU) ---

    // 1. Zjistíme, kolik dílů jsme v TOMTO konkrétním souboru reálně přidali
    int totalNewItems = 0;
    for (auto const& [id, aisObj] : m_display_map) {
        if (id >= fileId) totalNewItems++;
    }

    int currentItem = 0;

    for (auto const& [id, aisObj] : m_display_map) {
        // TOTO JE TA MAGIE: Přeskočíme staré díly z předchozích souborů
        if (id < fileId) continue;

        currentItem++;

        // Aktualizace UI pouze pro nové díly
        if (currentItem % 10 == 0 || currentItem == totalNewItems || currentItem == 1) {
            int percent = (currentItem * 100) / (totalNewItems > 0 ? totalNewItems : 1);
            updateStatus("Fáze 4/4: GPU a Výběr (" + std::to_string(percent) + " %) - Díl " + std::to_string(currentItem) + " z " + std::to_string(totalNewItems));
        }

        // 1. Zobrazení na obrazovce (jen nové díly)
        myContext->Display(aisObj, AIS_Shaded, 0, false);

        // 2. IHNED aktivujeme struktury pro výběr
        myContext->Activate(aisObj, 0);
        myContext->Activate(aisObj, AIS_Shape::SelectionMode(TopAbs_FACE));
    }

    // Staré generování výběru je smazáno (probíhá nyní ve Fázi 4)

    myContext->UpdateCurrentViewer();
    myView->FitAll();
    drawOriginAxes(); // Znovu nakreslí osy!
    updateStatus("Model úspěšně načten.");
    return true;
}

void OcctGtkGLAreaViewer::toggleClippingPlane()
{
    if (myView.IsNull() || myContext.IsNull()) return;

    // --- DETEKCE DVOJKLIKU POMOCÍ ČASOVAČE ---
    static gint64 last_click_time = 0;
    gint64 current_time = g_get_monotonic_time();
    // 400 000 mikrosekund (400 ms) je ideální doba pro zachycení běžného dvojkliku
    bool isDoubleClick = (current_time - last_click_time < 400000);
    last_click_time = isDoubleClick ? 0 : current_time; // Po dvojkliku časovač vyresetujeme

    // 1. Prvotní vytvoření roviny (pokud ještě neexistuje)
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

    // --- AKCE PRO DVOJKLIK (Otočení roviny) ---
    if (isDoubleClick) {
        gp_Pln currentPln = m_clipPlane->ToPlane();
        gp_Ax1 axis = currentPln.Axis();

        // Obrátíme normálu plochy
        axis.SetDirection(axis.Direction().Reversed());
        currentPln.SetAxis(axis);
        m_clipPlane->SetEquation(currentPln);

        m_clipPlane->SetOn(true); // Po otočení musí být řez vždy aktivní
        myView->Redraw();
        return;
    }

    // --- AKCE PRO JEDNOKLIK ---

    gp_Dir finalDir(0, 0, 1);
    gp_Pnt finalPos(0, 0, 0);
    bool hasExplicitDir = false;
    bool hasExplicitPos = false;
    gp_Dir fallbackDir(0, 0, 1);
    bool hasFallbackDir = false;

    // NOVÉ: Kontejner pro sběr vybraných bodů
    std::vector<gp_Pnt> selectedVertices;

    // 1. Zjištění výběru z plátna
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

            // NOVÉ: Sbíráme body do seznamu
            selectedVertices.push_back(pnt);

            finalPos = pnt;
            hasExplicitPos = true;
        }
    }

    // --- NOVÉ: Logika pro proložení roviny 3 body ---
    if (selectedVertices.size() == 3) {
        // Vytvoříme dva vektory z prvního bodu do zbylých dvou
        gp_Vec v1(selectedVertices[0], selectedVertices[1]);
        gp_Vec v2(selectedVertices[0], selectedVertices[2]);

        // Vektorový součin nám dá normálu (kolmici) k rovině
        gp_Vec normal = v1.Crossed(v2);

        // Kontrola, jestli body neleží v jedné přímce (normála by byla nulová)
        if (normal.Magnitude() > gp::Resolution()) {
            finalDir = gp_Dir(normal);
            finalPos = selectedVertices[0]; // Rovina prochází prvním bodem
            hasExplicitDir = true;
            hasExplicitPos = true;
            std::cout << "[ŘEZ] Vytvořena rovina ze 3 bodů." << std::endl;
        } else {
            std::cout << "[ŘEZ] CHYBA: Vybrané 3 body leží na jedné přímce!" << std::endl;
        }
    }
    // -----------------------------------------------
    else if (!hasExplicitDir && hasFallbackDir) {
        finalDir = fallbackDir;
        hasExplicitDir = true;
    }

    // 2. Vyhodnocení jednokliku
    if (hasExplicitDir || hasExplicitPos) {
        // A. Pokud má uživatel vybranou nějakou plochu/hranu, vytvoříme nový řez.
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
        // B. Pokud nemá nic vybráno, funguje tlačítko jako čistý Vypínač (Zapnout/Vypnout)
        bool isCurrentlyOn = m_clipPlane->IsOn();
        m_clipPlane->SetOn(!isCurrentlyOn);
    }

    myView->Redraw();
}

void OcctGtkGLAreaViewer::disableClippingPlane()
{
    // 1. Vypneme řez (pokud existuje)
    if (!m_clipPlane.IsNull()) {
        m_clipPlane->SetOn(false);
        m_clipState = 0;
    }

    // 2. Vyčistíme výběr
    if (!myContext.IsNull()) {
        myContext->ClearSelected(false);
    }

    // 3. Vynutíme překreslení okna
    if (!myView.IsNull()) {
        myView->Redraw();
    }
}

void OcctGtkGLAreaViewer::setShapeVisibility(int id, bool visible)
{
    if (myContext.IsNull()) return;

    bool screenNeedsUpdate = false;
    bool isMajorHide = false;

    // 1. SKRÝVÁNÍ CELÝCH HLAVNÍCH DÍLŮ
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
    // 2. SKRÝVÁNÍ POD-TVARŮ (MULTIBODY Tělesa i jednotlivé plošky)
    else if (m_subshape_map.find(id) != m_subshape_map.end()) {
        TopoDS_Shape subShape = m_subshape_map[id];

        // Hledání rodiče (hlavního grafického dílu) proplouváním vrstev
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

                // TADY JE TO KOUZLO:
                // Ať už klikneme na FACE nebo celé MULTIBODY těleso (SOLID),
                // rozložíme si náš požadavek na plošky a upravíme je všechny hromadně.
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

                // Překreslíme celou grafiku jen jednou!
                if (madeChanges) {
                    myContext->Redisplay(coloredObj, Standard_True);
                    screenNeedsUpdate = true;
                }

                // Pokud jsme skryli celé těleso (SOLID), zakážeme rekurzi do stromu,
                // abychom to nedělali 500x znovu pro každou dceřinou plochu
                if (subShape.ShapeType() == TopAbs_SOLID) {
                    isMajorHide = true;
                }
            }
        }
    }

    // 3. REKURZE DO STROMU
    if (m_tree_hierarchy.find(id) != m_tree_hierarchy.end()) {
        for (int childId : m_tree_hierarchy[id]) {
            if (isMajorHide && m_subshape_map.find(childId) != m_subshape_map.end()) {
                continue; // Tyto plošky jsme již zvládli vyřešit hromadně!
            }
            setShapeVisibility(childId, visible);
        }
    }

    // 4. PŘEKRESLENÍ MONITORU
    if (screenNeedsUpdate) {
        myContext->UpdateCurrentViewer();
        queue_draw();
    }
}

void OcctGtkGLAreaViewer::drawOriginAxes()
{
    if (myContext.IsNull()) return;

    gp_Pnt origin(0, 0, 0);

    // Vytvoření čar o délce 30mm
    Handle(AIS_Shape) aisX = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(30, 0, 0)));
    Handle(AIS_Shape) aisY = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(0, 30, 0)));
    Handle(AIS_Shape) aisZ = new AIS_Shape(BRepBuilderAPI_MakeEdge(origin, gp_Pnt(0, 0, 30)));

    // Zrušíme dynamické zvýrazňování (nechceme, aby osy šly vybírat myší)
    aisX->SetMutable(false);
    aisY->SetMutable(false);
    aisZ->SetMutable(false);

    // Nastavení barev
    aisX->SetColor(Quantity_NOC_RED);
    aisY->SetColor(Quantity_NOC_GREEN);
    aisZ->SetColor(Quantity_NOC_BLUE);

    // Nastavení tloušťky čar (zde tloušťka 2.0)
    aisX->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_RED, Aspect_TOL_SOLID, 2.0));
    aisY->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_GREEN, Aspect_TOL_SOLID, 2.0));
    aisZ->Attributes()->SetLineAspect(new Prs3d_LineAspect(Quantity_NOC_BLUE, Aspect_TOL_SOLID, 2.0));

    // Klíčový trik: Vrstva Topmost zajistí, že se osy vykreslí naposledy a prosvitnou přes každý model
    aisX->SetZLayer(Graphic3d_ZLayerId_Topmost);
    aisY->SetZLayer(Graphic3d_ZLayerId_Topmost);
    aisZ->SetZLayer(Graphic3d_ZLayerId_Topmost);

    // Zobrazíme na plátně (parametr -1 znamená, že nebudou aktivní pro výběr myší)
    myContext->Display(aisX, 0, -1, false);
    myContext->Display(aisY, 0, -1, false);
    myContext->Display(aisZ, 0, -1, false);
}

void OcctGtkGLAreaViewer::onRightClick(int n_press, double x, double y)
{
    // Zjistíme, jestli uživatel neklikl do prázdna
    int selectedId = getSelectedShapeId();

    if (selectedId != -1) {
        m_contextMenuTargetId = selectedId;

        // Nastavíme, aby menu vyskočilo přesně tam, kde je myš
        Gdk::Rectangle rect;
        rect.set_x(static_cast<int>(x));
        rect.set_y(static_cast<int>(y));
        rect.set_width(1);
        rect.set_height(1);

        m_contextMenu.set_pointing_to(rect);
        m_contextMenu.popup();
    }
}

