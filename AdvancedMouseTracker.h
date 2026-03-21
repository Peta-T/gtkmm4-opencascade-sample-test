#pragma once

#include <gtkmm.h>
#include <map>
#include <string>

class AdvancedMouseTracker
{
public:
    explicit AdvancedMouseTracker(Gtk::Widget* target_widget);
    virtual ~AdvancedMouseTracker();

    Gdk::ModifierType get_modifiers() const {
        return m_current_modifiers;
    }

protected:
    // Single Click
    virtual void on_MB1_cli() {}
    virtual void on_MB2_cli() {}
    virtual void on_MB3_cli() {}
    virtual void on_MB9_cli() {}
    virtual void on_MB8_cli() {}

    // Multi Click
    virtual void on_MB1_double_cli() {}
    virtual void on_MB1_triple_cli() {}
    virtual void on_MB1_quad_cli() {}
    virtual void on_MB2_double_cli() {}
    virtual void on_MB3_double_cli() {}
    virtual void on_MB9_double_cli() {}
    virtual void on_MB8_double_cli() {}

    // Release
    virtual void on_MB1_rel() {}
    virtual void on_MB2_rel() {}
    virtual void on_MB3_rel() {}
    virtual void on_MB9_rel() {}
    virtual void on_MB8_rel() {}

    // Drag Start/Update/End
    virtual void on_MB1_drag_start() {}
    virtual void on_MB1_drag_update(double dx, double dy) {}
    virtual void on_MB1_drag_end() {}

    virtual void on_MB2_drag_start() {}
    virtual void on_MB2_drag_update(double dx, double dy) {}
    virtual void on_MB2_drag_end() {}

    virtual void on_MB3_drag_start() {}
    virtual void on_MB3_drag_update(double dx, double dy) {}
    virtual void on_MB3_drag_end() {}

    virtual void on_MB9_drag_start() {}
    virtual void on_MB9_drag_update(double dx, double dy) {}
    virtual void on_MB9_drag_end() {}

    virtual void on_MB8_drag_start() {}
    virtual void on_MB8_drag_update(double dx, double dy) {}
    virtual void on_MB8_drag_end() {}

    // Combo Drag
    virtual void on_MB21_drag_start() {}
    virtual void on_MB21_drag_update(double dx, double dy) {}
    virtual void on_MB21_drag_end() {}

    virtual void on_MB23_drag_start() {}
    virtual void on_MB23_drag_update(double dx, double dy) {}
    virtual void on_MB23_drag_end() {}

    // Hover / Passive Motion
    virtual void on_mouse_move(double x, double y) {}

    // Scroll
    virtual void on_scroll(double dx, double dy) {}

protected:
    struct DragCandidate {
        bool is_pressed = false;
        bool valid_for_drag = true;
        double start_x = 0.0;
        double start_y = 0.0;
    };

    Gtk::Widget* m_target;
    const double DRAG_THRESHOLD = 5.0;

    Glib::RefPtr<Gtk::EventControllerMotion> m_controller_motion;
    Glib::RefPtr<Gtk::EventControllerLegacy> m_controller_legacy;
    Glib::RefPtr<Gtk::EventControllerScroll> m_controller_scroll;
    std::vector<Glib::RefPtr<Gtk::GestureClick>> m_gestures;

    std::map<guint, DragCandidate> m_drag_candidates;
    guint m_active_drag_id = 0;

    Gdk::ModifierType m_current_modifiers = (Gdk::ModifierType)0;

    void setup_gesture(guint btn_id);
    void on_mouse_motion_internal(double x, double y);
    void on_mouse_enter_internal(double x, double y);
    void on_mouse_leave_internal();
    bool on_scroll_internal(double dx, double dy);

    void on_mb_pressed_internal(int n_press, double x, double y, Gtk::GestureClick* gesture, guint button_id);
    void on_mb_released_internal(int n_press, double x, double y, Gtk::GestureClick* gesture, guint button_id);

    void on_mb_stopped_internal(guint button_id);
    bool on_event_legacy(const Glib::RefPtr<const Gdk::Event>& event);
    void force_drag_end(std::string reason);
};