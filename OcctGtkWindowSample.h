// File: OcctGtkWindowSample.h
#ifndef _OcctGtkWindowSample_HeaderFile
#define _OcctGtkWindowSample_HeaderFile

#include "OcctGtkGLAreaViewer.h"
#include <gtkmm.h>
#include <gtkmm/paned.h> // <--- Important: add this header

class OcctGtkWindowSample : public Gtk::Window
{
public:
  OcctGtkWindowSample();
  virtual ~OcctGtkWindowSample();

  void set_status(const Glib::ustring& message);
  void openFile(const Glib::ustring& filename);

protected:
  // Handlers
  void onFitAllClicked();
  void onViewFrontClicked();
  void onViewAxoClicked();

  void onOpenClicked();
  void on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog);

  void onTreeSelectionChanged();
  void onAboutClicked();

  // Keyboard handlers
  bool onModifiersChanged(Gdk::ModifierType ) { return myEventCtrlKey->forward(myViewer); }
  bool onKeyPressed(guint , guint , Gdk::ModifierType ) { return myEventCtrlKey->forward(myViewer); }
  void onKeyReleased(guint , guint , Gdk::ModifierType ) { myEventCtrlKey->forward(myViewer); }

protected:
  // Column definitions for TreeView
  class ModelColumns : public Gtk::TreeModel::ColumnRecord
  {
  public:
    ModelColumns()
    {
      add(m_col_name);
      add(m_col_type);
      add(m_col_id);
    }
    Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    Gtk::TreeModelColumn<Glib::ustring> m_col_type;
    Gtk::TreeModelColumn<int>           m_col_id;
  };

  ModelColumns m_columns;

protected:
  Gtk::Box    myVBox;

  // Toolbar
  Gtk::Box    myToolbar;
  Gtk::Button myBtnOpen;

  Gtk::Button myBtnFitAll;
  Gtk::Button myBtnViewFront;
  Gtk::Button myBtnViewAxo;

  Gtk::MenuButton myBtnMenu;
  Glib::RefPtr<Gio::SimpleActionGroup> m_action_group;
  Glib::RefPtr<Gio::Menu> m_menu;

  // Layout Containers
  Gtk::Box            m_content_box;

  // --- NEW: Paned (splitter) ---
  Gtk::Paned          m_paned;

  // TreeView
  Gtk::ScrolledWindow m_tree_scroller;
  Gtk::TreeView       m_tree_view;
  Glib::RefPtr<Gtk::TreeStore> m_tree_model;

  // Viewer
  OcctGtkGLAreaViewer myViewer;

  // Statusbar
  Gtk::Statusbar m_statusbar;

  Glib::RefPtr<Gtk::EventControllerKey> myEventCtrlKey;
};

#endif // _OcctGtkWindowSample_HeaderFile