// File: OcctGtkWindowSample.h
#ifndef _OcctGtkWindowSample_HeaderFile
#define _OcctGtkWindowSample_HeaderFile

#include "OcctGtkGLAreaViewer.h"
#include <gtkmm.h>
#include <gtkmm/paned.h> // Important for layout splitting

class OcctGtkWindowSample : public Gtk::Window
{
public:
  //! Main constructor
  OcctGtkWindowSample();
  virtual ~OcctGtkWindowSample();

  //! Updates the statusbar message
  void set_status(const Glib::ustring& message);

  //! Initiates loading of a STEP/XCAF file
  void openFile(const Glib::ustring& filename);

protected:
  // UI Event Handlers
  void onFitAllClicked();
  void onViewFrontClicked();
  void onViewAxoClicked();

  void onOpenClicked();
  void on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog);

  void onTreeSelectionChanged();
  void onToggleVisibility();
  void onAboutClicked();

  // Keyboard Event Handlers
  bool onModifiersChanged(Gdk::ModifierType) { return myEventCtrlKey->forward(myViewer); }
  bool onKeyPressed(guint theKeyVal, guint theKeyCode, Gdk::ModifierType theState);
  void onKeyReleased(guint, guint, Gdk::ModifierType) { myEventCtrlKey->forward(myViewer); }

  //! Focuses and selects a specific item in the browser tree by ID
  void focusTreeItemById(int searchId);

protected:
  // Column definitions for the Object Browser (TreeView)
  class ModelColumns : public Gtk::TreeModel::ColumnRecord
  {
  public:
    ModelColumns()
    {
      add(m_col_name);
      add(m_col_type);
      add(m_col_id);
      add(m_col_visible);
    }
    Gtk::TreeModelColumn<Glib::ustring> m_col_name;    // Object name
    Gtk::TreeModelColumn<Glib::ustring> m_col_type;    // Geometry type
    Gtk::TreeModelColumn<int>           m_col_id;      // Unique identifier
    Gtk::TreeModelColumn<bool>          m_col_visible; // Visibility state
  };

  ModelColumns m_columns;

protected:
  Gtk::Box    myVBox; // Main vertical container

  // Toolbar elements
  Gtk::Box    myToolbar;
  Gtk::Button myBtnOpen;

  Gtk::Button myBtnFitAll;
  Gtk::Button myBtnViewFront;
  Gtk::Button myBtnViewAxo;
  Gtk::Button m_btn_clip;

  Gtk::MenuButton myBtnMenu;
  Glib::RefPtr<Gio::SimpleActionGroup> m_action_group;
  Glib::RefPtr<Gio::Menu> m_menu;
  Gtk::PopoverMenu m_tree_popup;

  // Main Layout Containers
  Gtk::Box            m_content_box;
  Gtk::Paned          m_paned; // Splitter between tree and viewer

  // Object Browser (TreeView) elements
  Gtk::ScrolledWindow m_tree_scroller;
  Gtk::TreeView       m_tree_view;
  Glib::RefPtr<Gtk::TreeStore> m_tree_model;

  // 3D Viewport
  OcctGtkGLAreaViewer myViewer;

  // Application Status Bar
  Gtk::Statusbar m_statusbar;

  // Input Controllers
  Glib::RefPtr<Gtk::EventControllerKey> myEventCtrlKey;
};

#endif // _OcctGtkWindowSample_HeaderFile
