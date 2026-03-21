#include "OcctInputBridge.h"
#include "OcctGtkTools.h"

// --- CONDITIONAL COMPILATION FOR LOGGING ---
#ifdef OCCT_DEBUG_MOUSE
#include <iostream>
#include <iomanip>
#endif

// --- OCCT HEADERS ---
#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <StdSelect_ViewerSelector3d.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopAbs_ShapeEnum.hxx>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBndLib.hxx> // [NEW] For Bounding Box

// --- FOR PROPERTY CALCULATION ---
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <Extrema_ExtElC.hxx>
#include <Extrema_POnCurv.hxx>

#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Graphic3d_ZLayerId.hxx>

#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Sphere.hxx>
#include <gp_Cone.hxx>
#include <gp_Torus.hxx>
#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <Precision.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Bnd_Box.hxx> // [NEW]

#include <SelectMgr_EntityOwner.hxx>
#include <StdSelect_BRepOwner.hxx>

#include <sstream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =======================================================================
// AUXILIARY DEFINITIONS
// =======================================================================

enum SimplifiedType { GEOM_NONE, GEOM_POINT, GEOM_AXIS, GEOM_PLANE };

struct SimplifiedGeom {
    SimplifiedType type = GEOM_NONE;
    gp_Pnt point;
    gp_Lin axis;
    gp_Pln plane;
    TopoDS_Shape shape;
};

// =======================================================================
// AUXILIARY FUNCTIONS
// =======================================================================

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

static SimplifiedGeom analyzeShapeGeometry(const TopoDS_Shape& shape) {
    SimplifiedGeom res;
    res.shape = shape;

    if (shape.ShapeType() == TopAbs_VERTEX) {
        res.type = GEOM_POINT;
        res.point = BRep_Tool::Pnt(TopoDS::Vertex(shape));
    }
    else if (shape.ShapeType() == TopAbs_EDGE) {
        TopoDS_Edge edge = TopoDS::Edge(shape);
        BRepAdaptor_Curve curve(edge);
        if (curve.GetType() == GeomAbs_Line) {
            res.type = GEOM_AXIS;
            res.axis = curve.Line();
        }
        else if (curve.GetType() == GeomAbs_Circle) {
            res.type = GEOM_AXIS;
            res.axis = curve.Circle().Axis();
        }
    }
    else if (shape.ShapeType() == TopAbs_FACE) {
        TopoDS_Face face = TopoDS::Face(shape);
        BRepAdaptor_Surface surf(face);
        if (surf.GetType() == GeomAbs_Plane) {
            res.type = GEOM_PLANE;
            res.plane = surf.Plane();
        }
        else if (surf.GetType() == GeomAbs_Cylinder) {
            res.type = GEOM_AXIS;
            res.axis = surf.Cylinder().Axis();
        }
        else if (surf.GetType() == GeomAbs_Cone) {
            res.type = GEOM_AXIS;
            res.axis = surf.Cone().Axis();
        }
        else if (surf.GetType() == GeomAbs_Torus) {
            res.type = GEOM_AXIS;
            res.axis = surf.Torus().Axis();
        }
    }
    return res;
}

static std::string analyzeSingleShapeInfo(const TopoDS_Shape& shape, int id, const gp_Pnt& clickPnt)
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3);

    if (shape.IsNull()) return "";

    if (shape.ShapeType() == TopAbs_VERTEX) {
        gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(shape));
        ss << "POINT (" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
    }
    else if (shape.ShapeType() == TopAbs_EDGE) {
        TopoDS_Edge edge = TopoDS::Edge(shape);
        BRepAdaptor_Curve curve(edge);

        GProp_GProps LProps;
        BRepGProp::LinearProperties(edge, LProps);
        double length = LProps.Mass();

        double first = curve.FirstParameter();
        double last = curve.LastParameter();
        gp_Pnt pStart = curve.Value(first);
        gp_Pnt pEnd = curve.Value(last);

        if (curve.GetType() == GeomAbs_Line) {
            gp_Vec v(pStart, pEnd);
            ss << "LINE Len=" << length
               << " | dX=" << std::abs(v.X()) << " dY=" << std::abs(v.Y()) << " dZ=" << std::abs(v.Z());
        }
        else if (curve.GetType() == GeomAbs_Circle) {
            gp_Circ circ = curve.Circle();
            double r = circ.Radius();
            ss << (edge.Closed() ? "CIRCLE" : "ARC")
               << " D=" << (r * 2.0) << " R=" << r << " Len=" << length;
        }
        else {
            ss << "CURVE Len=" << length;
        }
    }
    else if (shape.ShapeType() == TopAbs_FACE) {
        TopoDS_Face face = TopoDS::Face(shape);
        BRepAdaptor_Surface surf(face);

        GProp_GProps SProps;
        BRepGProp::SurfaceProperties(shape, SProps);
        double area = SProps.Mass();

        if (surf.GetType() == GeomAbs_Plane) {
            gp_Pln plane = surf.Plane();
            gp_Dir norm = plane.Axis().Direction();
            double offset = plane.Distance(gp_Pnt(0,0,0));
            ss << "PLANAR FACE"
               << " | Normal(" << norm.X() << ", " << norm.Y() << ", " << norm.Z() << ")"
               << " | OriginDist=" << offset;
        }
        else if (surf.GetType() == GeomAbs_Cylinder) {
            gp_Cylinder cyl = surf.Cylinder();
            double r = cyl.Radius();
            gp_Dir axis = cyl.Axis().Direction();
            ss << "CYLINDER D=" << (r * 2.0) << " R=" << r
               << " | Axis(" << axis.X() << ", " << axis.Y() << ", " << axis.Z() << ")";
        }
        else if (surf.GetType() == GeomAbs_Sphere) {
            double r = surf.Sphere().Radius();
            ss << "SPHERE D=" << (r * 2.0) << " R=" << r;
        }
        else if (surf.GetType() == GeomAbs_Cone) {
            gp_Cone cone = surf.Cone();
            double semiAngleDeg = cone.SemiAngle() * (180.0 / M_PI);
            gp_Dir axis = cone.Axis().Direction();
            gp_Lin axisLine(cone.Axis());
            double rAtClick = axisLine.Distance(clickPnt);
            ss << "CONE D(at click)=" << (rAtClick * 2.0)
               << " Angle=" << semiAngleDeg << " deg"
               << " | Axis(" << axis.X() << ", " << axis.Y() << ", " << axis.Z() << ")";
        }
        else if (surf.GetType() == GeomAbs_Torus) {
            gp_Torus torus = surf.Torus();
            double majR = torus.MajorRadius();
            double minR = torus.MinorRadius();
            gp_Dir axis = torus.Axis().Direction();
            ss << "TORUS MajD=" << (majR * 2.0) << " MajR=" << majR
               << " | MinD=" << (minR * 2.0) << " MinR=" << minR
               << " | Axis(" << axis.X() << ", " << axis.Y() << ", " << axis.Z() << ")";
        }
        else {
            ss << "COMPLEX FACE";
        }
        ss << " | Area=" << area;
    }
    // --- SOLID / SHELL / COMPOUND ---
    else if (shape.ShapeType() == TopAbs_SOLID || shape.ShapeType() == TopAbs_SHELL || shape.ShapeType() == TopAbs_COMPSOLID || shape.ShapeType() == TopAbs_COMPOUND) {
        ss << shapeTypeToString(shape.ShapeType());

        GProp_GProps SProps;
        BRepGProp::SurfaceProperties(shape, SProps);
        double area = SProps.Mass();

        if (shape.ShapeType() == TopAbs_SOLID || shape.ShapeType() == TopAbs_COMPSOLID) {
            GProp_GProps VProps;
            BRepGProp::VolumeProperties(shape, VProps);
            double volume = VProps.Mass();

            // [NEW] Center of Gravity (COG)
            gp_Pnt cog = VProps.CentreOfMass();

            ss << " | Vol=" << volume << " | Mass(d=1)=" << volume;
            ss << " | COG(" << cog.X() << ", " << cog.Y() << ", " << cog.Z() << ")";
        }

        ss << " | Surf=" << area;

        // [NEW] Bounding Box (Envelope dimensions)
        Bnd_Box bbox;
        BRepBndLib::Add(shape, bbox);
        if (!bbox.IsVoid()) {
            double xmin, ymin, zmin, xmax, ymax, zmax;
            bbox.Get(xmin, ymin, zmin, xmax, ymax, zmax);
            ss << " | BBox(Dx=" << (xmax - xmin)
               << " Dy=" << (ymax - ymin)
               << " Dz=" << (zmax - zmin) << ")";
        }
    }
    else {
        ss << shapeTypeToString(shape.ShapeType());
    }

    ss << " | ID: " << id;
    return ss.str();
}

// =======================================================================
// CONSTRUCTOR
// =======================================================================

OcctInputBridge::OcctInputBridge(Gtk::Widget* target,
                                 const Handle(V3d_View)& view,
                                 const Handle(AIS_InteractiveContext)& context,
                                 AIS_ViewController* controller,
                                 std::function<void()> redraw_fn,
                                 std::function<void(const std::string&)> status_callback)
    : AdvancedMouseTracker(target),
      myView(view),
      myContext(context),
      myController(controller),
      myRedrawCallback(redraw_fn),
      myStatusCallback(status_callback),
      m_last_drag_x(0.0),
      m_last_drag_y(0.0),
      m_is_panning_active(false)
{
    m_controller_key = Gtk::EventControllerKey::create();
    m_controller_key->signal_key_pressed().connect(sigc::mem_fun(*this, &OcctInputBridge::on_key_pressed), false);
    m_controller_key->signal_key_released().connect(sigc::mem_fun(*this, &OcctInputBridge::on_key_released), false);
    target->add_controller(m_controller_key);
}

Aspect_VKeyFlags OcctInputBridge::get_occt_modifiers() {
    Aspect_VKeyFlags flags = OcctGtkTools::gtkModifiers2VKeys(get_modifiers());
    if (myController->Keys().IsKeyDown(Aspect_VKey_Shift))   flags |= Aspect_VKeyFlags_SHIFT;
    if (myController->Keys().IsKeyDown(Aspect_VKey_Control)) flags |= Aspect_VKeyFlags_CTRL;
    if (myController->Keys().IsKeyDown(Aspect_VKey_Alt))     flags |= Aspect_VKeyFlags_ALT;
    return flags;
}

Graphic3d_Vec2i OcctInputBridge::get_drag_pos(guint btn_id, double dx, double dy) {
    double abs_x = m_drag_candidates[btn_id].start_x + dx;
    double abs_y = m_drag_candidates[btn_id].start_y + dy;
    if (!myView->Window().IsNull()) {
        Graphic3d_Vec2d p(abs_x, abs_y);
        Graphic3d_Vec2d resultDouble = myView->Window()->ConvertPointToBacking(p) + Graphic3d_Vec2d(0.5);
        return Graphic3d_Vec2i((int)resultDouble.x(), (int)resultDouble.y());
    }
    return Graphic3d_Vec2i((int)abs_x, (int)abs_y);
}

// =======================================================================
// LINE VISUALIZATION
// =======================================================================

void OcctInputBridge::clearDistanceLine() {
    if (!m_dist_line.IsNull()) {
        myContext->Remove(m_dist_line, false);
        m_dist_line.Nullify();
    }
}

void OcctInputBridge::drawDistanceLine(const gp_Pnt& p1, const gp_Pnt& p2) {
    clearDistanceLine();

    if (p1.SquareDistance(p2) < Precision::SquareConfusion()) return;

    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(p1, p2);
    m_dist_line = new AIS_Shape(edge);

    m_dist_line->SetColor(Quantity_NOC_RED);
    m_dist_line->SetWidth(3.0);

    // Set TOPMOST layer (line visible over everything)
    m_dist_line->SetZLayer(Graphic3d_ZLayerId_Topmost);

    myContext->Display(m_dist_line, false);
    myContext->Deactivate(m_dist_line);
}

// =======================================================================
// MEASUREMENT
// =======================================================================
void OcctInputBridge::measureSelectedFeatures()
{
    if (myContext.IsNull()) return;

    clearDistanceLine();

    std::vector<TopoDS_Shape> selectedShapes;
    myContext->InitSelected();
    while (myContext->MoreSelected()) {
        Handle(SelectMgr_EntityOwner) owner = myContext->SelectedOwner();
        Handle(StdSelect_BRepOwner) brepOwner = Handle(StdSelect_BRepOwner)::DownCast(owner);
        if (!brepOwner.IsNull()) {
            selectedShapes.push_back(brepOwner->Shape());
        } else {
            Handle(AIS_InteractiveObject) obj = Handle(AIS_InteractiveObject)::DownCast(owner->Selectable());
            Handle(AIS_Shape) aisShape = Handle(AIS_Shape)::DownCast(obj);
            if (!aisShape.IsNull()) selectedShapes.push_back(aisShape->Shape());
        }
        myContext->NextSelected();
    }

    std::stringstream ss;
    ss << std::fixed << std::setprecision(3);

    if (selectedShapes.empty()) {
        if (myStatusCallback) myStatusCallback("Ready");
        return;
    }

    if (selectedShapes.size() == 1) {
        TopoDS_Shape shape = selectedShapes[0];
        int id = static_cast<int>(std::hash<TopoDS_Shape>{}(shape) & 0x7FFFFFFF);
        gp_Pnt p = myContext->MainSelector()->PickedPoint(1);
        ss << analyzeSingleShapeInfo(shape, id, p);
    }
    else if (selectedShapes.size() == 2) {
        const TopoDS_Shape& shape1 = selectedShapes[0];
        const TopoDS_Shape& shape2 = selectedShapes[1];

        BRepExtrema_DistShapeShape extrema(shape1, shape2);
        double minDist = -1.0;
        if (extrema.IsDone()) {
            minDist = extrema.Value();
            ss << "Min Dist: " << minDist;
        }

        SimplifiedGeom g1 = analyzeShapeGeometry(shape1);
        SimplifiedGeom g2 = analyzeShapeGeometry(shape2);

        // Auxiliary points for line drawing
        gp_Pnt pLine1, pLine2;
        bool hasLine = false;

        // --- 1. PLANE + POINT ---
        if ((g1.type == GEOM_PLANE && g2.type == GEOM_POINT) ||
            (g1.type == GEOM_POINT && g2.type == GEOM_PLANE))
        {
            const gp_Pln& pln = (g1.type == GEOM_PLANE) ? g1.plane : g2.plane;
            const gp_Pnt& pnt = (g1.type == GEOM_POINT) ? g1.point : g2.point;

            double u, v;
            ElSLib::Parameters(pln, pnt, u, v);
            pLine1 = pnt;
            pLine2 = ElSLib::Value(u, v, pln);
            hasLine = true;

            double perpDist = pln.Distance(pnt);
            ss << " | Perp Dist (Point-Plane): " << perpDist;
            if (perpDist < Precision::Confusion()) ss << " (Point is ON Plane)";
        }
        // --- 2. AXIS + POINT ---
        else if ((g1.type == GEOM_AXIS && g2.type == GEOM_POINT) ||
                 (g1.type == GEOM_POINT && g2.type == GEOM_AXIS))
        {
            const gp_Lin& lin = (g1.type == GEOM_AXIS) ? g1.axis : g2.axis;
            const gp_Pnt& pnt = (g1.type == GEOM_POINT) ? g1.point : g2.point;

            double distToAxis = lin.Distance(pnt);
            ss << " | Dist to Axis: " << distToAxis;

            double u = ElCLib::Parameter(lin, pnt);
            pLine1 = pnt;
            pLine2 = ElCLib::Value(u, lin);
            hasLine = true;
        }
        // --- 3. PLANE + PLANE ---
        else if (g1.type == GEOM_PLANE && g2.type == GEOM_PLANE) {
            double angleRad = g1.plane.Axis().Angle(g2.plane.Axis());
            double angleDeg = angleRad * (180.0 / M_PI);

            if (angleDeg < 0.1 || std::abs(angleDeg - 180.0) < 0.1) {
                double dist = g1.plane.Distance(g2.plane.Location());
                ss << " | Parallel Planes | Offset Dist: " << dist;

                pLine2 = g2.plane.Location();
                double u, v;
                ElSLib::Parameters(g1.plane, pLine2, u, v);
                pLine1 = ElSLib::Value(u, v, g1.plane);
                hasLine = true;
            } else {
                ss << " | Angle: " << angleDeg << " deg";
            }
        }
        // --- 4. AXIS + AXIS (SKEW / PARALLEL) ---
        else if (g1.type == GEOM_AXIS && g2.type == GEOM_AXIS) {
            double angleRad = g1.axis.Angle(g2.axis);
            double angleDeg = angleRad * (180.0 / M_PI);

            if (angleDeg < 0.1 || std::abs(angleDeg - 180.0) < 0.1) {
                ss << " | Parallel Axes | Dist: " << g1.axis.Distance(g2.axis);

                pLine2 = g2.axis.Location();
                double u = ElCLib::Parameter(g1.axis, pLine2);
                pLine1 = ElCLib::Value(u, g1.axis);
                hasLine = true;
            } else {
                ss << " | Angle: " << angleDeg << " deg";

                Extrema_ExtElC extremaAx(g1.axis, g2.axis, 1.0e-5);
                if (extremaAx.IsDone() && !extremaAx.IsParallel()) {
                    if (extremaAx.NbExt() > 0) {
                        double axisDist = std::sqrt(extremaAx.SquareDistance(1));
                        ss << " | Skew Axes Dist: " << axisDist;

                        Extrema_POnCurv P1_onCurve, P2_onCurve;
                        extremaAx.Points(1, P1_onCurve, P2_onCurve);
                        pLine1 = P1_onCurve.Value();
                        pLine2 = P2_onCurve.Value();
                        hasLine = true;
                    }
                }
            }
        }
        // --- 5. OTHERS (DEFAULT) ---
        else {
             // If we have min distance, use it as default
             if (minDist >= 0 && extrema.IsDone() && extrema.NbSolution() > 0) {
                 pLine1 = extrema.PointOnShape1(1);
                 pLine2 = extrema.PointOnShape2(1);
                 hasLine = true;
             }

             if ((g1.type == GEOM_PLANE && g2.type == GEOM_AXIS) ||
                 (g1.type == GEOM_AXIS && g2.type == GEOM_PLANE)) {
                const gp_Pln& pln = (g1.type == GEOM_PLANE) ? g1.plane : g2.plane;
                const gp_Lin& lin = (g1.type == GEOM_AXIS)  ? g1.axis  : g2.axis;
                double angleNormLine = pln.Axis().Direction().Angle(lin.Direction());
                double angleDeg = std::abs(90.0 - (angleNormLine * (180.0 / M_PI)));
                ss << " | Angle (Axis/Plane): " << angleDeg << " deg";
             }
        }

        // [NEW] Draw line and print dX/dY/dZ
        if (hasLine) {
            drawDistanceLine(pLine1, pLine2);

            // Delta values
            gp_Vec v(pLine1, pLine2);
            ss << " | dX=" << std::abs(v.X())
               << " dY=" << std::abs(v.Y())
               << " dZ=" << std::abs(v.Z());
        }
    }
    else {
        ss << "Selection: " << selectedShapes.size() << " items";
    }

    if (selectedShapes.size() >= 2) {
        double totalLen = 0.0;
        int countLinear = 0;
        for (const auto& s : selectedShapes) {
            if (s.ShapeType() == TopAbs_EDGE || s.ShapeType() == TopAbs_WIRE) {
                GProp_GProps LProps;
                BRepGProp::LinearProperties(s, LProps);
                totalLen += LProps.Mass();
                countLinear++;
            }
        }
        if (countLinear > 0) ss << " | Total Len=" << totalLen;
    }

    if (myStatusCallback) myStatusCallback(ss.str());
}

// =======================================================================
// EVENTS
// =======================================================================

void OcctInputBridge::on_mouse_move(double x, double y) {
    if (OcctGtkTools::gtkHandleMotionEvent(*myController, myView, Graphic3d_Vec2d(x, y), get_occt_modifiers())) {
        myRedrawCallback();
    }
}

void OcctInputBridge::on_scroll(double dx, double dy) {
    if (OcctGtkTools::gtkHandleScrollEvent(*myController, myView, Graphic3d_Vec2d(dx, dy))) {
        myRedrawCallback();
    }
}

void OcctInputBridge::on_MB1_cli() {
    double x = m_drag_candidates[GDK_BUTTON_PRIMARY].start_x;
    double y = m_drag_candidates[GDK_BUTTON_PRIMARY].start_y;

    Graphic3d_Vec2i pt;
    if (!myView->Window().IsNull()) {
         Graphic3d_Vec2d res = myView->Window()->ConvertPointToBacking(Graphic3d_Vec2d(x, y)) + Graphic3d_Vec2d(0.5);
         pt = Graphic3d_Vec2i((int)res.x(), (int)res.y());
    } else { pt = Graphic3d_Vec2i((int)x, (int)y); }

    myController->PressMouseButton(pt, Aspect_VKeyMouse_LeftButton, get_occt_modifiers(), false);
    myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_LeftButton, get_occt_modifiers(), false);

    if (!myContext.IsNull()) {
        Aspect_VKeyFlags mods = get_occt_modifiers();
        bool isMultiSelect = (mods & Aspect_VKeyFlags_SHIFT) || (mods & Aspect_VKeyFlags_CTRL);
        AIS_SelectionScheme scheme = isMultiSelect ? AIS_SelectionScheme_XOR : AIS_SelectionScheme_Replace;

        if (myContext->HasDetected()) {
            myContext->SelectDetected(scheme);
            myContext->UpdateCurrentViewer();
            measureSelectedFeatures();
        }
        else {
            if (!isMultiSelect) {
                myContext->ClearSelected(false);
                clearDistanceLine();
                if (myStatusCallback) myStatusCallback("Ready");
            }
        }
    }
    myRedrawCallback();
}

void OcctInputBridge::on_MB1_drag_start() {
    Graphic3d_Vec2i pt = get_drag_pos(GDK_BUTTON_PRIMARY, 0, 0);
    myController->PressMouseButton(pt, Aspect_VKeyMouse_LeftButton, get_occt_modifiers(), false);
}

void OcctInputBridge::on_MB1_drag_update(double dx, double dy) {
    Graphic3d_Vec2i pt = get_drag_pos(GDK_BUTTON_PRIMARY, dx, dy);
    myController->UpdateMousePosition(pt, Aspect_VKeyMouse_LeftButton, get_occt_modifiers(), false);
    myRedrawCallback();
}

void OcctInputBridge::on_MB1_drag_end() {
    Graphic3d_Vec2i pt = myController->LastMousePosition();
    myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_LeftButton, get_occt_modifiers(), false);
    myRedrawCallback();
}

void OcctInputBridge::on_MB2_drag_start() {
    m_is_panning_active = false;
    Graphic3d_Vec2i pt = get_drag_pos(GDK_BUTTON_MIDDLE, 0, 0);
   /// std::cout << "\n[LOG] 1. ROTATION START (MB2_drag_start) | Coordinates: " << pt.x() << ", " << pt.y() << std::endl;

    myController->PressMouseButton(pt, Aspect_VKeyMouse_MiddleButton, get_occt_modifiers(), false);
}

void OcctInputBridge::on_MB2_drag_update(double dx, double dy) {
    Graphic3d_Vec2i pt = get_drag_pos(GDK_BUTTON_MIDDLE, dx, dy);

    if (m_is_panning_active) {
     ///   std::cout << "[LOG] 4. RETURN TO ROTATION (MB2_drag_update after releasing MB3) | Coordinates: " << pt.x() << ", " << pt.y() << std::endl;

        myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_RightButton, get_occt_modifiers(), false);
        myController->PressMouseButton(pt, Aspect_VKeyMouse_MiddleButton, get_occt_modifiers(), false);
        m_is_panning_active = false;
    }

    myController->UpdateMousePosition(pt, Aspect_VKeyMouse_MiddleButton, get_occt_modifiers(), false);
    myRedrawCallback();
}

void OcctInputBridge::on_MB2_drag_end() {
    Graphic3d_Vec2i pt = myController->LastMousePosition();
   /// std::cout << "[LOG] 5. END OF ALL (MB2_drag_end) | Coordinates: " << pt.x() << ", " << pt.y() << std::endl;

    if (m_is_panning_active) {
        myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_RightButton, get_occt_modifiers(), false);
    } else {
        myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_MiddleButton, get_occt_modifiers(), false);
    }
    m_is_panning_active = false;
}

void OcctInputBridge::on_MB2_double_cli() {
    if (!myView.IsNull()) {
        myView->FitAll(0.01, false);
        myView->Invalidate();
        myRedrawCallback();
    }
}

void OcctInputBridge::on_MB23_drag_start() {
    // Not used directly in logic, monitored via update
}

void OcctInputBridge::on_MB23_drag_update(double dx, double dy) {
    Graphic3d_Vec2i pt = get_drag_pos(GDK_BUTTON_MIDDLE, dx, dy);

    if (!m_is_panning_active) {
       /// std::cout << "[LOG] 2. TRANSITION TO PANNING (MB23_drag_update - MB3 pressed) | Coordinates: " << pt.x() << ", " << pt.y() << std::endl;

        myController->ReleaseMouseButton(pt, Aspect_VKeyMouse_MiddleButton, get_occt_modifiers(), false);
        myController->PressMouseButton(pt, Aspect_VKeyMouse_RightButton, get_occt_modifiers(), false);
        m_is_panning_active = true;
    }

    // Print only occasionally to avoid flooding the console (every 15th pixel), but update constantly
    static int counter = 0;
    if (counter++ % 15 == 0) {
      ///  std::cout << "[LOG] 3. PANNING IN PROGRESS (MB23_drag_update) | Coordinates: " << pt.x() << ", " << pt.y() << std::endl;
    }

    myController->UpdateMousePosition(pt, Aspect_VKeyMouse_RightButton, get_occt_modifiers(), false);
    myRedrawCallback();
}

void OcctInputBridge::on_MB23_drag_end() {
   /// std::cout << "[LOG] MB23_drag_end (MB3 released - processed during next movement in MB2)" << std::endl;
}

void OcctInputBridge::on_MB3_cli() {}

bool OcctInputBridge::on_key_pressed(guint keyval, guint keycode, Gdk::ModifierType state) {
    Aspect_VKey vkey = OcctGtkTools::gtkKey2VKey(keyval, keycode);
    if (vkey != Aspect_VKey_UNKNOWN) {
        myController->KeyDown(vkey, myController->EventTime());
        myController->ProcessInput();
    }
    return false;
}

void OcctInputBridge::on_key_released(guint keyval, guint keycode, Gdk::ModifierType state) {
    Aspect_VKey vkey = OcctGtkTools::gtkKey2VKey(keyval, keycode);
    if (vkey != Aspect_VKey_UNKNOWN) {
        myController->KeyUp(vkey, myController->EventTime());
        myController->ProcessInput();
    }
}