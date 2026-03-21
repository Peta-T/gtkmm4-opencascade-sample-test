#ifndef _OcctInputBridge_HeaderFile
#define _OcctInputBridge_HeaderFile

#include "AdvancedMouseTracker.h"

#include <AIS_InteractiveContext.hxx>
#include <AIS_ViewController.hxx>
#include <V3d_View.hxx>
#include <AIS_Shape.hxx> // Required for m_dist_line

#include <functional>
#include <string>

// Forward declarations
class OcctInputBridge : public AdvancedMouseTracker
{
public:
    OcctInputBridge(Gtk::Widget* target,
                    const Handle(V3d_View)& view,
                    const Handle(AIS_InteractiveContext)& context,
                    AIS_ViewController* controller,
                    std::function<void()> redraw_fn,
                    std::function<void(const std::string&)> status_callback);

    virtual ~OcctInputBridge() = default;

protected:
    // Override methods from AdvancedMouseTracker
    virtual void on_MB1_cli() override;
    virtual void on_MB1_drag_start() override;
    virtual void on_MB1_drag_update(double dx, double dy) override;
    virtual void on_MB1_drag_end() override;

    virtual void on_MB2_drag_start() override;
    virtual void on_MB2_drag_update(double dx, double dy) override;
    virtual void on_MB2_drag_end() override;
    virtual void on_MB2_double_cli() override;

    virtual void on_MB23_drag_start() override;
    virtual void on_MB23_drag_update(double dx, double dy) override;
    virtual void on_MB23_drag_end() override;

    virtual void on_MB3_cli() override;

    virtual void on_mouse_move(double x, double y) override;
    virtual void on_scroll(double dx, double dy) override;

    // Keyboard
    bool on_key_pressed(guint keyval, guint keycode, Gdk::ModifierType state);
    void on_key_released(guint keyval, guint keycode, Gdk::ModifierType state);

private:
    Aspect_VKeyFlags get_occt_modifiers();
    Graphic3d_Vec2i get_drag_pos(guint btn_id, double dx, double dy);

    // Measurement logic
    void measureSelectedFeatures();

    // NEW: Line visualization
    void drawDistanceLine(const gp_Pnt& p1, const gp_Pnt& p2);
    void clearDistanceLine();

private:
    Handle(V3d_View) myView;
    Handle(AIS_InteractiveContext) myContext;
    AIS_ViewController* myController;
    std::function<void()> myRedrawCallback;
    std::function<void(const std::string&)> myStatusCallback;

    Glib::RefPtr<Gtk::EventControllerKey> m_controller_key;

    double m_last_drag_x;
    double m_last_drag_y;
    bool m_is_panning_active;

    // NEW: Storage for measurement visualization
    Handle(AIS_Shape) m_dist_line;
};

#endif // _OcctInputBridge_HeaderFile