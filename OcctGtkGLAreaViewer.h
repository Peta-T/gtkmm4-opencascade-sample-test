#ifndef _OcctGtkGLAreaViewer_HeaderFile
#define _OcctGtkGLAreaViewer_HeaderFile

#ifdef UNICODE
#undef UNICODE
#endif
#include <gtkmm.h>

#include <AIS_InteractiveContext.hxx>
#include <AIS_ViewController.hxx>
#include <AIS_ViewCube.hxx>
#include <AIS_Shape.hxx>
#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFApp_Application.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <XCAFPrs_AISObject.hxx>
#include <TDF_Label.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDataStd_Name.hxx>
#include <TopLoc_Location.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <SelectMgr_EntityOwner.hxx>
#include <StdSelect_BRepOwner.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <Graphic3d_ClipPlane.hxx>
#include <gtkmm/popovermenu.h>
#include <giomm/simpleactiongroup.h>
#include <gtkmm/gestureclick.h>

#include <map>
#include <string>
#include <vector>

class OcctInputBridge;

//! GTK GLArea widget with embedded OCCT Viewer.
class OcctGtkGLAreaViewer : public Gtk::GLArea, public AIS_ViewController
{
public:
  //! Main constructor.
  OcctGtkGLAreaViewer();
  virtual ~OcctGtkGLAreaViewer();

  const Handle(V3d_Viewer)& Viewer() const { return myViewer; }
  const Handle(V3d_View)& View() const { return myView; }
  const Handle(AIS_InteractiveContext)& Context() const { return myContext; }
  const TCollection_AsciiString& GetGlInfo() const { return myGlInfo; }
  void toggleClippingPlane();
  void disableClippingPlane();
  void drawOriginAxes();

  // New STEP/XCAF methods
  bool loadStepFile(const Glib::ustring& filePath,
                    Glib::RefPtr<Gtk::TreeStore>& treeModel,
                    const Gtk::TreeModelColumn<Glib::ustring>& colName,
                    const Gtk::TreeModelColumn<Glib::ustring>& colType,
                    const Gtk::TreeModelColumn<int>& colId,
                    const Gtk::TreeModelColumn<bool>& colVisible);

  void selectShapeById(int id);
  void fitToShapeById(int id);
  int getSelectedShapeId();
  void setShapeVisibility(int id, bool visible);

  // NEW: Signal for sending status messages to the window
  sigc::signal<void(const Glib::ustring&)> signal_status_message;
  sigc::signal<void(int)> signal_locate_in_tree;

  private:
    Handle(Graphic3d_ClipPlane) m_clipPlane;
    int m_clipState = 0;
    Gtk::PopoverMenu m_contextMenu;
    Glib::RefPtr<Gio::SimpleActionGroup> m_actionGroup;
    Glib::RefPtr<Gtk::GestureClick> m_rightClickGesture;
    int m_contextMenuTargetId = -1;

    // Zpracování pravého kliknutí
    void onRightClick(int n_press, double x, double y);

protected:
  // Input handling methods
  void updateModifiers();
  bool onModifiersChanged(Gdk::ModifierType theType);
  bool onKeyPressed(guint theKeyVal, guint theKeyCode, Gdk::ModifierType theType);
  void onKeyReleased(guint theKeyVal, guint theKeyCode, Gdk::ModifierType theType);
  void onMotionMove(double theX, double theY);
  void onMouseButtonPressed(int theNbPressed, double theX, double theY);
  void onMouseButtonReleased(int theNbPressed, double theX, double theY);
  bool onMouseScroll(double theDeltaX, double theDeltaY);

  void onGlAreaRealized();
  void onGlAreaReleased();
  bool onGlAreaRender(const Glib::RefPtr<Gdk::GLContext>& theGlCtx);

  void dumpGlInfo(bool theIsBasic, bool theToPrint);
  void initPixelScaleRatio();
  virtual void handleViewRedraw(const Handle(AIS_InteractiveContext)& theCtx,
                                const Handle(V3d_View)& theView) override;

  void activateSubShapeSelection();

protected:
  Handle(V3d_Viewer)             myViewer;
  Handle(V3d_View)               myView;
  Handle(AIS_InteractiveContext) myContext;
  Handle(AIS_ViewCube)           myViewCube;
  float                          myDevicePixelRatio = 1.0f;
  TCollection_AsciiString        myGlInfo;
  std::unique_ptr<OcctInputBridge> m_input_bridge;

  Glib::RefPtr<Gtk::EventControllerKey> myEventCtrlKey;
  Glib::RefPtr<Gtk::GestureClick>       myEventCtrlClick;
  Aspect_VKeyFlags myKeyModifiers = Aspect_VKeyFlags_NONE;

  Handle(TDocStd_Document) m_doc;
  Handle(AIS_Shape)        m_selected_part_ais;
  Handle(AIS_Shape)        m_ais_box;

  std::map<int, TDF_Label>       m_label_map;
  std::map<int, TopLoc_Location> m_loc_map;
  std::map<int, int>             m_depth_map;
  std::map<int, TopoDS_Shape>    m_subshape_map;
  std::map<int, Handle(AIS_InteractiveObject)> m_display_map;
  std::map<int, std::vector<int>> m_tree_hierarchy;

  int m_next_id = 1;
  bool m_is_wireframe = false;
};

#endif // _OcctGtkGLAreaViewer_HeaderFile
