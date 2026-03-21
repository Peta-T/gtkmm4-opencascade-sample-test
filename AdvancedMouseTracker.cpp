#include "AdvancedMouseTracker.h"
#include <cmath>
#include <iostream>
#include <gdk/gdk.h>

namespace {
    // Helper function for masks
    Gdk::ModifierType get_mask_for_btn(guint btn_id) {
        switch (btn_id) {
            case GDK_BUTTON_PRIMARY: return Gdk::ModifierType::BUTTON1_MASK;
            case GDK_BUTTON_MIDDLE: return Gdk::ModifierType::BUTTON2_MASK;
            case GDK_BUTTON_SECONDARY: return Gdk::ModifierType::BUTTON3_MASK;
            default: return static_cast<Gdk::ModifierType>(0);
        }
    }
}

AdvancedMouseTracker::AdvancedMouseTracker(Gtk::Widget* target_widget)
    : m_target(target_widget)
{
    if (!m_target) return;

    // 1. Motion
    m_controller_motion = Gtk::EventControllerMotion::create();
    m_controller_motion->set_propagation_phase(Gtk::PropagationPhase::CAPTURE);
    m_controller_motion->signal_motion().connect(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mouse_motion_internal));
    m_controller_motion->signal_enter().connect(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mouse_enter_internal));
    m_controller_motion->signal_leave().connect(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mouse_leave_internal));
    m_target->add_controller(m_controller_motion);

    // 2. Scroll
    m_controller_scroll = Gtk::EventControllerScroll::create();
    m_controller_scroll->set_flags(Gtk::EventControllerScroll::Flags::VERTICAL);
    m_controller_scroll->signal_scroll().connect(sigc::mem_fun(*this, &AdvancedMouseTracker::on_scroll_internal), true);
    m_target->add_controller(m_controller_scroll);

    // 3. Gestures
    setup_gesture(GDK_BUTTON_PRIMARY);
    setup_gesture(GDK_BUTTON_MIDDLE);
    setup_gesture(GDK_BUTTON_SECONDARY);
    setup_gesture(8);
    setup_gesture(9);

    // 4. Legacy (Fix MB8/MB9)
    m_controller_legacy = Gtk::EventControllerLegacy::create();
    m_controller_legacy->set_propagation_phase(Gtk::PropagationPhase::CAPTURE);
    m_controller_legacy->signal_event().connect(sigc::mem_fun(*this, &AdvancedMouseTracker::on_event_legacy), false);
    m_target->add_controller(m_controller_legacy);
}

AdvancedMouseTracker::~AdvancedMouseTracker() {}

void AdvancedMouseTracker::setup_gesture(guint btn_id) {
    auto g = Gtk::GestureClick::create();
    g->set_button(btn_id);
    g->set_propagation_phase(Gtk::PropagationPhase::CAPTURE);

    g->signal_pressed().connect(
        sigc::bind(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mb_pressed_internal), g.get(), btn_id)
    );
    g->signal_released().connect(
        sigc::bind(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mb_released_internal), g.get(), btn_id)
    );

    g->signal_stopped().connect(sigc::bind(sigc::mem_fun(*this, &AdvancedMouseTracker::on_mb_stopped_internal), btn_id));

    m_target->add_controller(g);
    m_gestures.push_back(g);
    m_drag_candidates[btn_id] = DragCandidate{};
}

void AdvancedMouseTracker::on_mb_pressed_internal(int n_press, double x, double y, Gtk::GestureClick* gesture, guint button_id)
{
    if (gesture) {
        m_current_modifiers = gesture->get_current_event_state();
    }

    if (n_press == 1) {
        if (m_drag_candidates[button_id].is_pressed) return;
        m_drag_candidates[button_id].is_pressed = true;
        m_drag_candidates[button_id].valid_for_drag = true;
        m_drag_candidates[button_id].start_x = x;
        m_drag_candidates[button_id].start_y = y;
    }

    switch (button_id) {
        case GDK_BUTTON_PRIMARY:
            if (n_press == 2) on_MB1_double_cli();
            else if (n_press == 3) on_MB1_triple_cli();
            else if (n_press == 4) on_MB1_quad_cli();
            break;
        case GDK_BUTTON_MIDDLE:
            if (n_press == 2) on_MB2_double_cli();
            break;
        case GDK_BUTTON_SECONDARY:
            if (n_press == 2) on_MB3_double_cli();
            break;
        case 9:
            if (n_press == 2) on_MB9_double_cli();
            break;
        case 8:
            if (n_press == 2) on_MB8_double_cli();
            break;
    }
}

void AdvancedMouseTracker::on_mouse_motion_internal(double x, double y)
{
    m_current_modifiers = m_controller_motion->get_current_event_state();
    Gdk::ModifierType state = m_controller_motion->get_current_event_state();
    for (auto& [btn_id, candidate] : m_drag_candidates) {
        if (candidate.is_pressed) {
            Gdk::ModifierType mask = get_mask_for_btn(btn_id);
            if ((int)mask != 0 && (state & mask) != mask) {
                candidate.is_pressed = false;
                if (m_active_drag_id != 0) force_drag_end("Safety Check Fail");
                return;
            }
        }
    }

    if (m_active_drag_id != 0)
    {
        auto& c = m_drag_candidates[m_active_drag_id];
        double dx = x - c.start_x;
        double dy = y - c.start_y;

        if (m_active_drag_id == GDK_BUTTON_PRIMARY)          on_MB1_drag_update(dx, dy);
        else if (m_active_drag_id == GDK_BUTTON_SECONDARY)  on_MB3_drag_update(dx, dy);
        else if (m_active_drag_id == 9)                     on_MB9_drag_update(dx, dy);
        else if (m_active_drag_id == 8)                     on_MB8_drag_update(dx, dy);
        else if (m_active_drag_id == GDK_BUTTON_MIDDLE) {
             bool left_down = m_drag_candidates[GDK_BUTTON_PRIMARY].is_pressed;
             bool right_down = m_drag_candidates[GDK_BUTTON_SECONDARY].is_pressed;
             if (right_down)       on_MB23_drag_update(dx, dy);
             else if (left_down)  on_MB21_drag_update(dx, dy);
             else                 on_MB2_drag_update(dx, dy);
        }
    }
    else
    {
        bool drag_started_now = false;
        for (auto& [btn_id, candidate] : m_drag_candidates) {
            if (candidate.is_pressed && candidate.valid_for_drag) {
                double dist = std::sqrt(std::pow(x - candidate.start_x, 2) + std::pow(y - candidate.start_y, 2));

                if (dist > DRAG_THRESHOLD) {
                    m_active_drag_id = btn_id;
                    drag_started_now = true;

                    if (btn_id == GDK_BUTTON_PRIMARY)           on_MB1_drag_start();
                    else if (btn_id == GDK_BUTTON_SECONDARY)    on_MB3_drag_start();
                    else if (btn_id == 9)                       on_MB9_drag_start();
                    else if (btn_id == 8)                       on_MB8_drag_start();
                    else if (btn_id == GDK_BUTTON_MIDDLE) {
                          bool left_down = m_drag_candidates[GDK_BUTTON_PRIMARY].is_pressed;
                          bool right_down = m_drag_candidates[GDK_BUTTON_SECONDARY].is_pressed;
                          if (right_down)       on_MB23_drag_start();
                          else if (left_down)  on_MB21_drag_start();
                          else                 on_MB2_drag_start();
                    }
                    break;
                }
            }
        }

        if (!drag_started_now) {
            on_mouse_move(x, y);
        }
    }
}

void AdvancedMouseTracker::on_mb_released_internal(int n_press, double x, double y, Gtk::GestureClick* gesture, guint button_id)
{
    if (gesture) {
        m_current_modifiers = gesture->get_current_event_state();
    }

    bool was_dragging_owner = (m_active_drag_id == button_id);
    bool is_any_drag_active = (m_active_drag_id != 0);

    if (n_press == 1 && !is_any_drag_active) {
        switch (button_id) {
            case GDK_BUTTON_PRIMARY:   on_MB1_cli(); break;
            case GDK_BUTTON_MIDDLE:    on_MB2_cli(); break;
            case GDK_BUTTON_SECONDARY: on_MB3_cli(); break;
            case 9:                    on_MB9_cli(); break;
            case 8:                    on_MB8_cli(); break;
        }
    }

    switch (button_id) {
        case GDK_BUTTON_PRIMARY:   on_MB1_rel(); break;
        case GDK_BUTTON_MIDDLE:    on_MB2_rel(); break;
        case GDK_BUTTON_SECONDARY: on_MB3_rel(); break;
        case 9:                    on_MB9_rel(); break;
        case 8:                    on_MB8_rel(); break;
    }

    if (was_dragging_owner) {
        if (button_id == GDK_BUTTON_PRIMARY)          on_MB1_drag_end();
        else if (button_id == GDK_BUTTON_SECONDARY)  on_MB3_drag_end();
        else if (button_id == 9)                      on_MB9_drag_end();
        else if (button_id == 8)                      on_MB8_drag_end();
        else if (button_id == GDK_BUTTON_MIDDLE)      on_MB2_drag_end();

        m_active_drag_id = 0;
    }

    m_drag_candidates[button_id].is_pressed = false;
    m_drag_candidates[button_id].valid_for_drag = true;
}

void AdvancedMouseTracker::on_mb_stopped_internal(guint button_id)
{
    if ((button_id == 9 || button_id == 8) && m_drag_candidates[button_id].is_pressed) {
        if (m_active_drag_id == button_id) return;
        return;
    }
}

void AdvancedMouseTracker::force_drag_end(std::string reason)
{
    if (m_active_drag_id != 0) {
        if (m_active_drag_id == 1) on_MB1_drag_end();
        else if (m_active_drag_id == 2) on_MB2_drag_end();
        else if (m_active_drag_id == 3) on_MB3_drag_end();
        else if (m_active_drag_id == 9) on_MB9_drag_end();
        else if (m_active_drag_id == 8) on_MB8_drag_end();
        m_active_drag_id = 0;
    }
    for (auto& [id, c] : m_drag_candidates) {
        if (c.is_pressed) c.valid_for_drag = false;
    }
}

bool AdvancedMouseTracker::on_event_legacy(const Glib::RefPtr<const Gdk::Event>& event)
{
    GdkEvent* c_event = const_cast<GdkEvent*>(event->gobj());
    if (gdk_event_get_event_type(c_event) == GDK_BUTTON_RELEASE)
    {
        guint btn_id = gdk_button_event_get_button(c_event);
        if ((btn_id == 8 || btn_id == 9) && m_drag_candidates.count(btn_id) && m_drag_candidates[btn_id].is_pressed)
        {
            on_mb_released_internal(1, 0, 0, nullptr, btn_id);
            return false;
        }
    }
    return false;
}

bool AdvancedMouseTracker::on_scroll_internal(double dx, double dy) {
    m_current_modifiers = m_controller_scroll->get_current_event_state();
    on_scroll(dx, dy);
    return true;
}

void AdvancedMouseTracker::on_mouse_enter_internal(double, double) {}
void AdvancedMouseTracker::on_mouse_leave_internal() {
    if (m_active_drag_id != 0) force_drag_end("Mouse left window");
}