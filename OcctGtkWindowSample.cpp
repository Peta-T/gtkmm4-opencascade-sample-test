// File: OcctGtkWindowSample.cpp
#include "OcctGtkWindowSample.h"

#include <AIS_Shape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Message.hxx>
#include <Standard_Version.hxx>
#include <iostream>

#include <Graphic3d_CLight.hxx>
#include <gp_Dir.hxx>
#include <Quantity_Color.hxx>
#include "app_icon.h"
#include <gdk/gdkkeysyms.h>

// ================================================================
// Function : OcctGtkWindowSample
// ================================================================
OcctGtkWindowSample::OcctGtkWindowSample()
: myVBox(Gtk::Orientation::VERTICAL),
  myToolbar(Gtk::Orientation::HORIZONTAL)
{
  set_title("gstepview");
  set_default_size(1024, 768);

  myVBox.set_spacing(0);
  set_child(myVBox);

  // 1. TOOLBAR SETUP
  myToolbar.set_spacing(5);
  myToolbar.set_margin(5);

  // OPEN Button
  myBtnOpen.set_icon_name("document-open");
  myBtnOpen.set_tooltip_text("Open STEP File");
  myBtnOpen.signal_clicked().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onOpenClicked));
  myToolbar.append(myBtnOpen);

  // Fit All Button
  myBtnFitAll.set_icon_name("zoom-fit-best");
  myBtnFitAll.set_tooltip_text("Fit All");
  myBtnFitAll.signal_clicked().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onFitAllClicked));
  myToolbar.append(myBtnFitAll);

  // Front Button
  myBtnViewFront.set_label("Front");
  myBtnViewFront.signal_clicked().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onViewFrontClicked));
  myToolbar.append(myBtnViewFront);

  // Iso Button
  myBtnViewAxo.set_label("Iso");
  myBtnViewAxo.signal_clicked().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onViewAxoClicked));
  myToolbar.append(myBtnViewAxo);

  // Toggle Clip Plane Button
  m_btn_clip.set_label("Toggle Clip Plane");
  m_btn_clip.signal_clicked().connect([this]() {
      myViewer.toggleClippingPlane();
  });
  myToolbar.append(m_btn_clip);

  // 1. Empty flexible box that pushes the menu completely to the right
  Gtk::Box* spacer = Gtk::manage(new Gtk::Box(Gtk::Orientation::HORIZONTAL));
  spacer->set_hexpand(true);
  myToolbar.append(*spacer);

  // 2. Action registration for "About" click
  m_action_group = Gio::SimpleActionGroup::create();
  m_action_group->add_action("about", sigc::mem_fun(*this, &OcctGtkWindowSample::onAboutClicked));
  insert_action_group("win", m_action_group); // "win" prefix for window actions

  // 3. Creation of the dropdown menu itself
  m_menu = Gio::Menu::create();
  m_menu->append("About", "win.about");

  // 4. Button setup (3 lines = "open-menu-symbolic")
  myBtnMenu.set_icon_name("open-menu-symbolic");
  myBtnMenu.set_menu_model(m_menu);
  myBtnMenu.set_tooltip_text("Application Menu");
  myToolbar.append(myBtnMenu);

  myVBox.append(myToolbar);
  myVBox.append(*Gtk::manage(new Gtk::Separator(Gtk::Orientation::HORIZONTAL)));

  // 2. CONTENT AREA (Paned Layout)
  m_content_box.set_orientation(Gtk::Orientation::HORIZONTAL);
  m_content_box.set_expand(true);

  m_paned.set_orientation(Gtk::Orientation::HORIZONTAL);
  m_paned.set_expand(true);
  m_paned.set_position(250);

  // -- TreeView --
  m_tree_model = Gtk::TreeStore::create(m_columns);
  m_tree_view.set_model(m_tree_model);

  int cols_count = m_tree_view.append_column("Object Browser", m_columns.m_col_name);
  if(auto* pCol = m_tree_view.get_column(cols_count - 1)) {
      pCol->set_expand(true);
      // Link the "sensitive" property (graying out) with our visibility column
      auto cell = dynamic_cast<Gtk::CellRendererText*>(pCol->get_first_cell());
      if (cell) pCol->add_attribute(cell->property_sensitive(), m_columns.m_col_visible);
  }


  m_tree_view.append_column("Type", m_columns.m_col_type);
  m_tree_view.set_enable_tree_lines(true);
  m_tree_view.signal_cursor_changed().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onTreeSelectionChanged));

  // Dummy data (initial root)
  auto row = *(m_tree_model->append());
  row[m_columns.m_col_name] = "Scene Root";
  row[m_columns.m_col_type] = "Group";
  row[m_columns.m_col_id]   = 0;
  row[m_columns.m_col_visible] = true;

  m_tree_scroller.set_child(m_tree_view);
  m_tree_scroller.set_policy(Gtk::PolicyType::AUTOMATIC, Gtk::PolicyType::AUTOMATIC);

  // -- Viewer --
  myViewer.set_hexpand(true);
  myViewer.set_vexpand(true);

  // Insert into Paned
  m_paned.set_start_child(m_tree_scroller);
  m_paned.set_end_child(myViewer);

  m_content_box.append(m_paned);
  myVBox.append(m_content_box);

  // 3. STATUSBAR
  myVBox.append(*Gtk::manage(new Gtk::Separator(Gtk::Orientation::HORIZONTAL)));
  m_statusbar.set_margin(5);
  myVBox.append(m_statusbar);

  set_status("Ready");

  // Controller setup
  myEventCtrlKey = Gtk::EventControllerKey::create();
  myEventCtrlKey->signal_modifiers().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onModifiersChanged), false);
  myEventCtrlKey->signal_key_pressed().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onKeyPressed), false);
  myEventCtrlKey->signal_key_released().connect(sigc::mem_fun(*this, &OcctGtkWindowSample::onKeyReleased), false);
  add_controller(myEventCtrlKey);

  // Initialize Viewer Axes
  {
    myViewer.drawOriginAxes();
  }

  if (!myViewer.View().IsNull())
  {
      Handle(V3d_Viewer) aViewer = myViewer.View()->Viewer();

      // 1. Remove all previous lights
      while (!aViewer->DefinedLights().IsEmpty())
      {
          aViewer->DelLight(aViewer->DefinedLights().First());
      }
      myViewer.View()->SetLightOff();

      // 2. AMBIENT LIGHT (Base CAD look)
      // Ensures that shadows are not pure black, but gray.
      Handle(Graphic3d_CLight) aAmb = new Graphic3d_CLight(Graphic3d_TOLS_AMBIENT);
      aAmb->SetColor(Quantity_NOC_WHITE);
      aAmb->SetIntensity(0.40f); // 40% light is omnipresent
      aViewer->AddLight(aAmb);
      myViewer.View()->SetLightOn(aAmb);

      // 3. WEAK HEADLIGHT (For gloss/highlights)
      // Very weak light from the camera, so surfaces perpendicular to the view
      // get a little shiny, but not "overexposed".
      Handle(Graphic3d_CLight) aFill = new Graphic3d_CLight(Graphic3d_TOLS_DIRECTIONAL);
      aFill->SetDirection(gp_Dir(-0.30, -0.30, -1.0)); // From the camera perspective
      aFill->SetColor(Quantity_NOC_WHITE);
      aFill->SetIntensity(0.020f); // Only 2%
      aFill->SetHeadlight(true);   // Moves with the camera
      aViewer->AddLight(aFill);
      myViewer.View()->SetLightOn(aFill);

      // Update view
      myViewer.View()->Update();
  }

  myViewer.signal_status_message.connect(sigc::mem_fun(*this, &OcctGtkWindowSample::set_status));
  // When the canvas sends a locate signal, call our method
  myViewer.signal_locate_in_tree.connect(sigc::mem_fun(*this, &OcctGtkWindowSample::focusTreeItemById));

  // Tree Context Menu setup
  auto tree_action_group = Gio::SimpleActionGroup::create();
  tree_action_group->add_action("toggle_visibility", sigc::mem_fun(*this, &OcctGtkWindowSample::onToggleVisibility));
  m_tree_view.insert_action_group("tree", tree_action_group);

  auto tree_menu = Gio::Menu::create();
  tree_menu->append("Hide / Show", "tree.toggle_visibility");

  m_tree_popup.set_parent(m_tree_view);
  m_tree_popup.set_menu_model(tree_menu);
  m_tree_popup.set_has_arrow(false);

  // Capture Right Click on TreeView
  auto click_controller = Gtk::GestureClick::create();
  click_controller->set_button(GDK_BUTTON_SECONDARY);
  click_controller->signal_pressed().connect([this](int n_press, double x, double y) {

      Gtk::TreeModel::Path path;
      Gtk::TreeViewColumn* column = nullptr;
      int cell_x, cell_y;

      // 1. Check what is exactly under the mouse
      bool on_item = m_tree_view.get_path_at_pos((int)x, (int)y, path, column, cell_x, cell_y);
      if (on_item) {
          // If GTK hit the item, select it immediately
          m_tree_view.get_selection()->select(path);
      }

      // 2. Ultimate fallback: Check what is currently selected in the tree
      auto iter = m_tree_view.get_selection()->get_selected();

      if (iter) {
          path = m_tree_model->get_path(iter);
          column = m_tree_view.get_column(0);

          if (column) {
              Gdk::Rectangle cell_rect;
              m_tree_view.get_cell_area(path, *column, cell_rect);
              // Direct menu exactly to the vertical center of the row
              int safe_y = cell_rect.get_y() + (cell_rect.get_height() / 2);
              m_tree_popup.set_pointing_to(Gdk::Rectangle((int)x, safe_y, 1, 1));
          } else {
              m_tree_popup.set_pointing_to(Gdk::Rectangle((int)x, (int)y, 1, 1));
          }

          m_tree_popup.set_position(Gtk::PositionType::RIGHT);
          m_tree_popup.popup();
      }
  });
  m_tree_view.add_controller(click_controller);
}

OcctGtkWindowSample::~OcctGtkWindowSample()
{
}

bool OcctGtkWindowSample::onKeyPressed(guint theKeyVal, guint theKeyCode, Gdk::ModifierType theState)
{
    // Global emergency brake for ESC
    if (theKeyVal == GDK_KEY_Escape) {
        myViewer.disableClippingPlane();
    }

    // Forward the signal using the built-in GTK mechanism
    return myEventCtrlKey->forward(myViewer);
}

void OcctGtkWindowSample::openFile(const Glib::ustring& filename)
{
  set_status("Loading: " + filename);

  bool res = myViewer.loadStepFile(filename, m_tree_model, m_columns.m_col_name, m_columns.m_col_type, m_columns.m_col_id, m_columns.m_col_visible);

  if (res) {
      set_status("Loaded: " + filename);
  } else {
      set_status("Failed to load: " + filename);
  }
}

void OcctGtkWindowSample::set_status(const Glib::ustring& message)
{
  m_statusbar.remove_all_messages();
  m_statusbar.push(message);
}

// --- VIEW HANDLERS ---
void OcctGtkWindowSample::onFitAllClicked()
{
  if (!myViewer.View().IsNull())
  {
    myViewer.View()->FitAll(0.01, false);
    myViewer.View()->Invalidate();
    myViewer.queue_draw();
    set_status("View: Fit All");
  }
}

void OcctGtkWindowSample::onViewFrontClicked()
{
  if (!myViewer.View().IsNull())
  {
    myViewer.View()->SetProj(V3d_Yneg);
    myViewer.View()->FitAll(0.01, false);
    myViewer.View()->Invalidate();
    myViewer.queue_draw();
    set_status("View: Front Projection");
  }
}

void OcctGtkWindowSample::onViewAxoClicked()
{
  if (!myViewer.View().IsNull())
  {
    myViewer.View()->SetProj(V3d_XposYnegZpos);
    myViewer.View()->FitAll(0.01, false);
    myViewer.View()->Invalidate();
    myViewer.queue_draw();
    set_status("View: Axonometric Projection");
  }
}

// --- FILE HANDLERS ---

void OcctGtkWindowSample::onOpenClicked()
{
  auto dialog = Gtk::FileDialog::create();
  dialog->set_title("Open STEP File");

  auto filters = Gio::ListStore<Gtk::FileFilter>::create();
  auto filterStep = Gtk::FileFilter::create();
  filterStep->set_name("STEP Files");
  filterStep->add_pattern("*.step");
  filterStep->add_pattern("*.stp");
  filters->append(filterStep);
  dialog->set_filters(filters);

  dialog->open(sigc::bind(sigc::mem_fun(*this, &OcctGtkWindowSample::on_file_dialog_finish), dialog));
}

void OcctGtkWindowSample::on_file_dialog_finish(const Glib::RefPtr<Gio::AsyncResult>& result, const Glib::RefPtr<Gtk::FileDialog>& dialog)
{
  try
  {
    auto file = dialog->open_finish(result);
    openFile(file->get_path());
  }
  catch (const Gtk::DialogError& err) {
    set_status("Open cancelled.");
  }
  catch (const Glib::Error& err) {
    set_status("Error: " + Glib::ustring(err.what()));
  }
}

void OcctGtkWindowSample::onTreeSelectionChanged()
{
  auto selection = m_tree_view.get_selection();
  auto iter = selection->get_selected();
  if(iter)
  {
    Glib::ustring name = (*iter)[m_columns.m_col_name];
    int id = (*iter)[m_columns.m_col_id];

    // Call the viewer method for highlighting
    myViewer.selectShapeById(id);

    set_status("Selected object: " + name + " (ID: " + std::to_string(id) + ")");
  }
}

void OcctGtkWindowSample::onAboutClicked(){
  Gtk::AboutDialog* dialog = new Gtk::AboutDialog();
  dialog->set_transient_for(*this);
  dialog->set_modal(true);

  dialog->set_program_name("gstepview OCCT Gtkmm 3D Viewer");
  dialog->set_version("1.0");
  dialog->set_copyright("© 2026");

  dialog->set_comments("Application for viewing and measuring STEP models.\n\n"
                       "Rotate model by middle mouse button\n\n"
                       "Pan model by middle and right mouse button\n\n"
                       "Zoom by mouse wheel rotate\n\n\n\n"
                       "Author: PetaT petaemail@seznam.cz\n\n"
                       "Special thanks to Mr. Kirill Gavrilov "
                       "for providing the original source codes:\n"
                       "• OcctGlTools.h / OcctGlTools.cpp\n"
                       "• OcctGtkTools.h / OcctGtkTools.cpp");

  try {
      auto bytes = Glib::Bytes::create(app_icon_png, app_icon_png_len);
      auto logo = Gdk::Texture::create_from_bytes(bytes);
      dialog->set_logo(logo);
  } catch (...) {
      dialog->set_logo_icon_name("applications-graphics");
  }

  dialog->signal_hide().connect([dialog]() { delete dialog; });
  dialog->present();
}

void OcctGtkWindowSample::onToggleVisibility()
{
    auto iter = m_tree_view.get_selection()->get_selected();
    if (iter) {
        bool current_vis = (*iter)[m_columns.m_col_visible];
        bool new_vis = !current_vis;

        std::function<void(Gtk::TreeRow&)> update_tree =
            [&](Gtk::TreeRow& row) {
                row[m_columns.m_col_visible] = new_vis;
                for (auto child : row.children()) {
                    update_tree(child);
                }
            };

        update_tree(*iter);

        // Perform hide/show in OpenCASCADE memory
        int id = (*iter)[m_columns.m_col_id];
        myViewer.setShapeVisibility(id, new_vis);

        // Recalculate viewer once after all changes are done
        if (!myViewer.Context().IsNull()) {
            myViewer.Context()->UpdateCurrentViewer();
        }
        myViewer.queue_draw();

        set_status(new_vis ? "Object is now visible" : "Object hidden");
    }
}

void OcctGtkWindowSample::focusTreeItemById(int searchId)
{
    if (!m_tree_model) return;

    std::function<bool(Gtk::TreeModel::Children)> searchLevel = [&](Gtk::TreeModel::Children children) {
        for (auto iter = children.begin(); iter != children.end(); ++iter) {

            Gtk::TreeModel::Row row = *iter;
            int rowId = row[m_columns.m_col_id];

            if (rowId == searchId) {
                // FOUND IT!
                Gtk::TreePath path = m_tree_model->get_path(iter);

                m_tree_view.expand_to_path(path);
                m_tree_view.get_selection()->select(path);
                m_tree_view.scroll_to_row(path);

                set_status("Found in tree: " + Glib::ustring(row[m_columns.m_col_name]));
                return true;
            }

            if (searchLevel(row.children())) return true;
        }
        return false;
    };

    // Start searching from the root
    searchLevel(m_tree_model->children());
}
