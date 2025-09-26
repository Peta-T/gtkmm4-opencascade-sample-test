/* build :
g++ -std=c++17 -g main.cc -o app $(pkg-config --cflags --libs gtkmm-4.0 epoxy librsvg-2.0) -I/usr/include/gtk-4.0 -I/usr/include/opencascade -I/usr/include/opencascade/Standard -I/mingw64/include/gtk-4.0 -I/mingw64/include/opencascade -I/mingw64/include/opencascade/Standard -I/usr/include/gtk-4.0/gdk/x11 -L/usr/local/lib/x86_64-linux-gnu -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep -lTKCAF -lTKCDF -lTKDCAF -lTKDECascade -lTKDEGLTF -lTKDEIGES -lTKDEOBJ -lTKDEPLY -lTKDE -lTKDESTEP -lTKDESTL -lTKDEVRML -lTKDraw -lTKernel -lTKExpress -lTKFeat -lTKFillet -lTKG2d -lTKG3d -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKLCAF -lTKMath -lTKMesh -lTKMeshVS -lTKOffset -lTKOpenGl -lTKOpenGlTest -lTKPrim -lTKQADraw -lTKRWMesh -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKTObjDRAW -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest -lTKXCAF -lTKXDEDRAW -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj -lTKXmlXCAF -lTKXSBase -lTKXSDRAWDE -lTKXSDRAWGLTF -lTKXSDRAWIGES -lTKXSDRAWOBJ -lTKXSDRAWPLY -lTKXSDRAW -lTKXSDRAWSTEP -lTKXSDRAWSTL -lTKXSDRAWVRML  -D_USE_MATH_DEFINES
*/

#include <stdio.h>



#include <gtkmm.h>
#include <epoxy/gl.h>

//#undef GLAPI // Undefine GLAPI to avoid redefinition issues with epoxy

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp> // For glm::rotate
#include <glm/gtx/transform.hpp>    // For glm::translate
#include <vector>
#include <string>
#include <filesystem>
#include <map>
#include <list> // For std::list in STEPImporter
#include <limits> // For std::numeric_limits
#include <glm/gtc/quaternion.hpp> // For glm::quat
#include <glm/gtx/quaternion.hpp> // For glm::toMat4, glm::slerp, glm::angleAxis

// OpenCASCADE Includes for BRep and Meshing
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx> // Added for processWire
#include <TopoDS_Compound.hxx> // Added for processComp
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColStd_Array1OfInteger.hxx> // Corrected include directive
#include <gp_Pnt.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Geom_Curve.hxx>
#include <TopLoc_Location.hxx> // Added for transformations
#include <gp_Trsf.hxx>         // Added for transformations
#include <gp_XYZ.hxx>          // Added for transformations
#include <Precision.hxx>       // Added for Precision::Confusion
#include <Quantity_Color.hxx>  // Added for colors
#include <BRepTools_WireExplorer.hxx> // Added for wire processing
#include <ShapeAnalysis_Edge.hxx>     // Added for edge analysis
#include <BRepAdaptor_Curve.hxx>      // Added for curve adaptation
#include <gp_Circ.hxx>                // Added for circle type in curves
#include <GCPnts_TangentialDeflection.hxx> // Added for edge discretization
#include <Poly.hxx>                   // Added for Poly::ComputeNormals

// XCAF for STEP import
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <Standard_Version.hxx>
#include <IFSelect_ReturnStatus.hxx> // Added for STEP import status
#include <Interface_Static.hxx>      // Added for STEP import settings
#include <STEPCAFControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_Label.hxx>       // Added for XCAF labels
#include <TDF_LabelSequence.hxx> // Added for XCAF label sequences
#include <TDF_ChildIterator.hxx> // Added for XCAF label iteration

// OpenCASCADE Geometry Types for surface analysis
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAbs_SurfaceType.hxx> // Required for GeomAbs_SurfaceType enum values

// GDK specific includes for native window handle (platform-dependent)
#ifdef GDK_WINDOWING_X11
#include <gdk/x11/gdkx.h>
#endif
//#ifdef GDK_WINDOWING_WIN32
//#include <gdk/win32/gdkwin32.h>
//#endif

// Custom comparator for gp_Pnt to use it in std::map
struct PointComparator {
    bool operator()(const gp_Pnt& p1, const gp_Pnt& p2) const {
        if (p1.X() != p2.X()) return p1.X() < p2.X();
        if (p1.Y() != p2.Y()) return p1.Y() < p2.Y();
        return p1.Z() < p2.Z();
    }
};

struct Color {
    float r, g, b, a;
    Color(float _r, float _g = 0.0f, float _b = 0.0f, float _a = 1.0f) : r(_r), g(_g), b(_b), a(_a) {}
};

const char* vertex_shader_src = R"(
#version 130
in vec4 in_position;
in vec4 in_color;
uniform mat4 mvp;
flat out vec4 v_color;

void main() {
    gl_Position = mvp * in_position;
    v_color = in_color;
}
)";

const char* fragment_shader_src = R"(
#version 130
flat in vec4 v_color;
out vec4 frag_color;
uniform bool is_line_pass;
uniform vec4 line_color_uniform;
uniform bool apply_normal_correction_tint; // New uniform
uniform vec4 normal_correction_tint_color; // New uniform

void main() {
    vec4 final_color;
    if (is_line_pass) {
        final_color = line_color_uniform; // For lines (black) and highlights (orange/magenta)
    } else {
        // For faces
        if (apply_normal_correction_tint) {
            // Mix original vertex color with the tint color
            final_color = mix(v_color, normal_correction_tint_color, 0.4); // Increased tint intensity to 40%
        } else {
            final_color = v_color;
        }
    }
    frag_color = final_color;
}
)";

// Helper function to check OpenGL errors
void GL_CHECK_ERROR_FUNC(const char* file, int line) {
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << "OpenGL error in " << file << ":" << line << " - Code: " << err << std::endl;
    }
}
#define GL_CHECK_ERROR GL_CHECK_ERROR_FUNC(__FILE__, __LINE__)

// --- Dune3D Structures and STEP Importer (from provided example) ---
namespace Dune3D {

// util/fs_util.hpp content
std::string path_to_string(const std::filesystem::path& p) {
    return p.string();
}

// shapes.hpp content
struct Vertex {
    float x, y, z;
    // Add a constructor to initialize members
    Vertex(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) : x(_x), y(_y), z(_z) {}

    bool operator<(const Vertex& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
    Vertex operator+(const Vertex& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Vertex operator/(float scalar) const { return {x / scalar, y / scalar, z / scalar}; }
    void operator+=(const Vertex& other) { x += other.x; y += other.y; z += other.z; }
    void operator/=(float scalar) { x /= scalar; y /= scalar; z /= scalar; }
};

struct FaceColor {
    float r, g, b;
};

struct Face {
    std::vector<Vertex> vertices;
    std::vector<Vertex> normals; // Not used in current shader, but good to keep
    std::vector<glm::ivec3> triangle_indices; // Indices into 'vertices' for this face
    FaceColor color;
    GeomAbs_SurfaceType surface_type = GeomAbs_OtherSurface; // New: Store the surface type
};

using Edge = std::vector<Vertex>; // A sequence of vertices forming an edge

struct Result {
    std::vector<Face> faces;
    std::vector<Edge> edges;
    // No need for a custom destructor if members are std::vector
};

struct Shapes {
    std::vector<TopoDS_Shape> shapes;
};

namespace STEPImporter {

#define USER_PREC (0.14)
#define USER_ANGLE (0.52359878) // M_PI / 6

#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
#define HORIZON_NEW_OCC
#endif

class STEPImporter {
public:
    STEPImporter(const std::filesystem::path& filename);
    bool is_loaded() const { return loaded; }
    Result get_faces_and_points();
    std::vector<TopoDS_Shape> get_shapes();

private:
    bool readSTEP(const char* fname);
    void processWire(const TopoDS_Wire& wire, const glm::dmat4& mat);
    bool processFace(const TopoDS_Face& face, Quantity_Color* color, const glm::dmat4& mat);
    bool processShell(const TopoDS_Shape& shape, Quantity_Color* color, const glm::dmat4& mat);
    bool getColor(TDF_Label label, Quantity_Color& color);
    bool processSolid(const TopoDS_Shape& shape, const glm::dmat4& mat_in = glm::dmat4(1.0));
    bool processComp(const TopoDS_Shape& shape, const glm::dmat4& mat_in = glm::dmat4(1.0));
    bool processNode(const TopoDS_Shape& shape);

    Handle(XCAFApp_Application) m_app;
    Handle(TDocStd_Document) m_doc;
    Handle(XCAFDoc_ShapeTool) m_assy;
    Handle(XCAFDoc_ColorTool) m_color;
    bool loaded = false;
    bool hasSolid = false; // Used internally by processSolid

    Result* result = nullptr; // Pointer to the result object being filled
};

STEPImporter::STEPImporter(const std::filesystem::path& filename)
{
    m_app = XCAFApp_Application::GetApplication();

    m_app->NewDocument("MDTV-XCAF", m_doc);
    if (!readSTEP(Dune3D::path_to_string(filename).c_str())) {
        std::cerr << "error loading " << filename << std::endl;
        loaded = false;
        return;
    }
    std::cerr << "loaded " << filename << std::endl;
    loaded = true;

    m_assy = XCAFDoc_DocumentTool::ShapeTool(m_doc->Main());
    m_color = XCAFDoc_DocumentTool::ColorTool(m_doc->Main());
}

bool STEPImporter::readSTEP(const char *fname)
{
    STEPCAFControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(fname);

    if (stat != IFSelect_RetDone)
        return false;

    // Enable user-defined shape precision
    if (!Interface_Static::SetIVal("read.precision.mode", 1))
        return false;

    // Set the shape conversion precision to USER_PREC (default 0.0001 has too
    // many triangles)
    if (!Interface_Static::SetRVal("read.precision.val", USER_PREC))
        return false;

    // set other translation options
    reader.SetColorMode(true);  // use model colors
    reader.SetNameMode(false);  // don't use label names
    reader.SetLayerMode(false); // ignore LAYER data

    if (!reader.Transfer(m_doc)) {
        m_doc->Close();
        return false;
    }

    // are there any shapes to translate?
    if (reader.NbRootsForTransfer() < 1)
        return false;

    return true;
}

void STEPImporter::processWire(const TopoDS_Wire &wire, const glm::dmat4 &mat)
{
    // This function was commented out in the previous version as it was not directly
    // used for rendering edges in the main drawing loop, and its original purpose
    // of populating `result->points` was unclear given the `Result` struct definition.
    // Edge data is now generated in `get_faces_and_points` using `GCPnts_TangentialDeflection`.
    // Keeping it commented out.
    /*
    for (BRepTools_WireExplorer expl(wire); expl.More(); expl.Next()) {
        const auto &edge = expl.Current();
        ShapeAnalysis_Edge sh_an;
        auto first = sh_an.FirstVertex(edge);
        {
            gp_Pnt pnt = BRep_Tool::Pnt(first);
            glm::dvec4 gpnt(pnt.X(), pnt.Y(), pnt.Z(), 1);
            auto pt = mat * gpnt;
            // result->points.emplace_back(pt.x, pt.y, pt.z);
        }

        auto curve = BRepAdaptor_Curve(edge);
        const auto curvetype = curve.GetType();
        if (curvetype == GeomAbs_Circle) {
            gp_Circ c = curve.Circle();
            auto pnt = c.Location();
            {
                glm::dvec4 gpnt(pnt.X(), pnt.Y(), pnt.Z(), 1);
                auto pt = mat * gpnt;
                // result->points.emplace_back(pt.x, pt.y, pt.z);
            }
        }
    }
    */
}

bool STEPImporter::processFace(const TopoDS_Face &face, Quantity_Color *color, const glm::dmat4 &mat)
{
    if (Standard_True == face.IsNull())
        return false;

    TopLoc_Location loc;
    Standard_Boolean isTessellate(Standard_False);
    Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);

    if (triangulation.IsNull() || triangulation->Deflection() > USER_PREC + Precision::Confusion())
        isTessellate = Standard_True;

    if (isTessellate) {
        BRepMesh_IncrementalMesh IM(face, USER_PREC, Standard_False, USER_ANGLE);
        triangulation = BRep_Tool::Triangulation(face, loc);
    }

    if (triangulation.IsNull() == Standard_True)
        return false;

    Quantity_Color lcolor;

    // check for a face color; this has precedence over SOLID colors
    do {
        TDF_Label L;

        if (m_color->ShapeTool()->Search(face, L)) {
            if (m_color->GetColor(L, XCAFDoc_ColorGen, lcolor) || m_color->GetColor(L, XCAFDoc_ColorCurv, lcolor)
                || m_color->GetColor(L, XCAFDoc_ColorSurf, lcolor))
                color = &lcolor;
        }
    } while (0);

    Poly::ComputeNormals(triangulation);

#ifndef HORIZON_NEW_OCC
    const TColgp_Array1OfPnt &arrPolyNodes = triangulation->Nodes();
    const Poly_Array1OfTriangle &arrTriangles = triangulation->Triangulation();
    const TShort_Array1OfShortReal &arrNormals = triangulation->Normals();
#endif


    result->faces.emplace_back();
    auto &face_out = result->faces.back();

    // Determine and store the surface type
    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    if (!surface.IsNull()) {
        if (surface->IsKind(STANDARD_TYPE(Geom_Plane))) {
            face_out.surface_type = GeomAbs_Plane;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_CylindricalSurface))) {
            face_out.surface_type = GeomAbs_Cylinder;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_ConicalSurface))) {
            face_out.surface_type = GeomAbs_Cone;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_SphericalSurface))) {
            face_out.surface_type = GeomAbs_Sphere;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_ToroidalSurface))) {
            face_out.surface_type = GeomAbs_Torus;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_BezierSurface))) {
            face_out.surface_type = GeomAbs_BezierSurface;
        } else if (surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface))) {
            face_out.surface_type = GeomAbs_BSplineSurface;
        } else {
            // For other complex surfaces like Geom_SurfaceOfRevolution or Geom_SurfaceOfLinearExtrusion,
            // which don't have a direct GeomAbs_ enum, default to OtherSurface.
            face_out.surface_type = GeomAbs_OtherSurface;
        }
    } else {
        face_out.surface_type = GeomAbs_OtherSurface;
    }


    if (color) {
        face_out.color.r = color->Red();
        face_out.color.g = color->Green();
        face_out.color.b = color->Blue();
    }
    else {
        face_out.color = {0.0f, 1.0f, 0.0f}; // Default bright green if no color
    }
    face_out.vertices.reserve(triangulation->NbNodes());

    std::map<Vertex, std::vector<size_t>> pts_map;
    for (int i = 1; i <= triangulation->NbNodes(); i++) {
#ifdef HORIZON_NEW_OCC
        gp_XYZ v(triangulation->Node(i).Coord());
#else
        gp_XYZ v(arrPolyNodes(i).Coord());
#endif
        const glm::dvec4 vg(v.X(), v.Y(), v.Z(), 1.0);
        const auto vt = mat * vg;
        // Use the new constructor for Vertex
        const Dune3D::Vertex vertex(static_cast<float>(vt.x), static_cast<float>(vt.y), static_cast<float>(vt.z));
        pts_map[vertex].push_back(i - 1);
        face_out.vertices.push_back(vertex);
    }

    face_out.normals.reserve(triangulation->NbNodes());
    for (int i = 1; i <= triangulation->NbNodes(); i++) {
#ifdef HORIZON_NEW_OCC
        const auto n = triangulation->Normal(i);
        glm::dvec4 vg(n.X(), n.Y(), n.Z(), 0.0);
#else
        auto offset = (i - 1) * 3 + 1;
        auto x = arrNormals(offset + 0);
        auto y = arrNormals(offset + 1);
        auto z = arrNormals(offset + 2);
        glm::dvec4 vg(x, y, z, 0.0);
#endif
        auto vt = mat * vg;
        vt = glm::normalize(vt); // Normalize after transformation
        // Use the new constructor for Vertex
        face_out.normals.emplace_back(static_cast<float>(vt.x), static_cast<float>(vt.y), static_cast<float>(vt.z));
    }

    // average normals at coincident vertices
    for (const auto &[k, v] : pts_map) {
        if (v.size() > 1) {
            Dune3D::Vertex n_acc(0.0f, 0.0f, 0.0f); // Initialize with floats
            for (const auto idx : v) {
                n_acc += face_out.normals.at(idx);
            }
            n_acc /= static_cast<float>(v.size());
            for (const auto idx : v) {
                face_out.normals.at(idx) = n_acc;
            }
        }
    }

    face_out.triangle_indices.reserve(triangulation->NbTriangles());
    for (int i = 1; i <= triangulation->NbTriangles(); i++) {
        int a, b, c;
#ifdef HORIZON_NEW_OCC
        triangulation->Triangle(i).Get(a, b, c);
#else
        arrTriangles(i).Get(a, b, c);
#endif
        // Check face orientation and reverse winding order if necessary
        if (face.Orientation() == TopAbs_FORWARD) {
            face_out.triangle_indices.emplace_back(a - 1, b - 1, c - 1);
        } else { // TopAbs_REVERSED
            face_out.triangle_indices.emplace_back(a - 1, c - 1, b - 1);
        }
    }

    return true;
}

bool STEPImporter::processShell(const TopoDS_Shape &shape, Quantity_Color *color, const glm::dmat4 &mat)
{
    TopoDS_Iterator it;
    bool ret = false;

    for (it.Initialize(shape, false, false); it.More(); it.Next()) {
        const TopoDS_Face &face = TopoDS::Face(it.Value());

        if (processFace(face, color, mat))
            ret = true;
    }

    return ret;
}

bool STEPImporter::getColor(TDF_Label label, Quantity_Color &color)
{
    while (true) {
        if (m_color->GetColor(label, XCAFDoc_ColorGen, color))
            return true;
        else if (m_color->GetColor(label, XCAFDoc_ColorSurf, color))
            return true;
        else if (m_color->GetColor(label, XCAFDoc_ColorCurv, color))
            return true;

        label = label.Father();

        if (label.IsNull())
            break;
    };

    return false;
}

bool STEPImporter::processSolid(const TopoDS_Shape &shape, const glm::dmat4 &mat_in)
{
    TDF_Label label = m_assy->FindShape(shape, Standard_False);
    bool ret = false;

    hasSolid = true;

    Quantity_Color col;
    Quantity_Color *lcolor = NULL;

    if (label.IsNull()) {
    }
    else {
        if (getColor(label, col))
            lcolor = &col;
    }

    TopLoc_Location loc = shape.Location();
    gp_Trsf T = loc.Transformation();
    gp_XYZ coord = T.TranslationPart();

    auto mat = mat_in * glm::translate(glm::dvec3(coord.X(), coord.Y(), coord.Z()));

    gp_XYZ axis;
    Standard_Real angle;

    if (T.GetRotation(axis, angle)) {
        glm::dvec3 gaxis(axis.X(), axis.Y(), axis.Z());
        double angle_f = angle;
        mat = glm::rotate(mat, angle_f, gaxis);
    }

    TopoDS_Iterator it;
    for (it.Initialize(shape, false, false); it.More(); it.Next()) {
        const TopoDS_Shape &subShape = it.Value();

        if (processShell(subShape, lcolor, mat))
            ret = true;
    }

    return ret;
}


bool STEPImporter::processComp(const TopoDS_Shape &shape, const glm::dmat4 &mat_in)
{
    TopoDS_Iterator it;
    TopLoc_Location loc = shape.Location();
    gp_Trsf T = loc.Transformation();
    gp_XYZ coord = T.TranslationPart();

    auto mat = mat_in * glm::translate(glm::dvec3(coord.X(), coord.Y(), coord.Z()));

    gp_XYZ axis;
    Standard_Real angle;

    if (T.GetRotation(axis, angle)) {
        glm::dvec3 gaxis(axis.X(), axis.Y(), axis.Z());
        double angle_f = angle;
        mat = glm::rotate(mat, angle_f, gaxis);
    }

    bool ret = false;

    for (it.Initialize(shape, false, false); it.More(); it.Next()) {
        const TopoDS_Shape &subShape = it.Value();
        TopAbs_ShapeEnum stype = subShape.ShapeType();

        hasSolid = false;

        switch (stype) {
        case TopAbs_COMPOUND:
        case TopAbs_COMPSOLID:
            if (processComp(subShape, mat))
                ret = true;
            break;

        case TopAbs_SOLID:
            if (processSolid(subShape, mat))
                ret = true;
            break;

        case TopAbs_SHELL:
            if (processShell(subShape, NULL, mat))
                ret = true;
            break;

        case TopAbs_FACE:
            if (processFace(TopoDS::Face(subShape), NULL, mat))
                ret = true;
            break;

        default:
            break;
        }
    }

    return ret;
}

bool STEPImporter::processNode(const TopoDS_Shape &shape)
{
    TopAbs_ShapeEnum stype = shape.ShapeType();
    bool ret = true;
    // Define an identity matrix for initial calls where no specific transformation is needed
    glm::dmat4 identity_mat(1.0);

    switch (stype) {
    case TopAbs_COMPOUND:
    case TopAbs_COMPSOLID:
        if (processComp(shape, identity_mat)) // Pass identity_mat
            ret = true;
        break;

    case TopAbs_SOLID:
        if (processSolid(shape, identity_mat)) // Pass identity_mat
            ret = true;
        break;

    case TopAbs_SHELL:
        if (processShell(shape, NULL, identity_mat)) // Pass identity_mat
            ret = true;
        break;

    case TopAbs_FACE:
        if (processFace(TopoDS::Face(shape), NULL, identity_mat)) // Pass identity_mat
            ret = true;
        break;

    default:
            break;
        }
    return ret;
}

Result STEPImporter::get_faces_and_points()
{
    Result res;
    result = &res; // Set the member pointer to the local result

    TDF_LabelSequence frshapes;
    m_assy->GetFreeShapes(frshapes);

    int nshapes = frshapes.Length();
    int id = 1;
    std::cerr << "shapes " << nshapes << std::endl;
    while (id <= nshapes) {
        TopoDS_Shape shape = m_assy->GetShape(frshapes.Value(id));
        if (!shape.IsNull() && processNode(shape)) {
            TopExp_Explorer topex(shape, TopAbs_EDGE);
            std::list<TopoDS_Shape> edges_processed; // Keep track of processed edges to avoid duplicates
            while (topex.More()) {
                auto edge = TopoDS::Edge(topex.Current());
                bool skip = false;
                for (const auto &other : edges_processed) {
                    if (other.IsSame(topex.Current())) {
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    topex.Next();
                    continue;
                }
                edges_processed.push_back(edge); // Mark as processed
                {
                    auto curve = BRepAdaptor_Curve(edge);
                    // M_PI / 16 is approx 11.25 degrees, 1e3 is max points
                    GCPnts_TangentialDeflection discretizer(curve, M_PI / 16.0, 1e3); // Use a reasonable deflection and max points, ensure float literal
                    auto &e = res.edges.emplace_back(); // Add a new Edge (vector of Vertices)
                    if (discretizer.NbPoints() > 0) {
                        int nbPoints = discretizer.NbPoints();
                        for (int i = 1; i <= nbPoints; i++) {
                            const gp_Pnt pnt = discretizer.Value(i);
                            e.emplace_back(static_cast<float>(pnt.X()), static_cast<float>(pnt.Y()), static_cast<float>(pnt.Z()));
                        }
                    }
                }
                topex.Next();
            }
        }
        ++id;
    }

    result = nullptr; // Clear the member pointer
    return res;
}

std::vector<TopoDS_Shape> STEPImporter::get_shapes()
{
    std::vector<TopoDS_Shape> r;
    TDF_LabelSequence frshapes;
    m_assy->GetFreeShapes(frshapes);

    int nshapes = frshapes.Length();
    int id = 1;
    std::cerr << "shapes " << nshapes << std::endl;
    while (id <= nshapes) {
        TopoDS_Shape shape = m_assy->GetShape(frshapes.Value(id));
        if (!shape.IsNull()) {
            r.push_back(shape);
        }
        ++id;
    }
    return r;
}

Result import(const std::filesystem::path& filename)
{
    STEPImporter importer(filename);
    if (!importer.is_loaded())
        return {};
    return importer.get_faces_and_points();
}

std::shared_ptr<const Shapes> import_shapes(const std::filesystem::path& filename)
{
    STEPImporter importer(filename);
    if (!importer.is_loaded())
        return {};
    auto shapes = std::make_shared<Shapes>();
    shapes->shapes = importer.get_shapes();
    return shapes;
}

} // namespace STEPImporter
} // namespace Dune3D

// New structures for highlighting
enum HighlightType { NONE, FACE, EDGE, VERTEX };

struct HighlightedElement {
    HighlightType type = NONE;
    int id = -1; // Index of the highlighted face/edge/vertex
};

class MyGLArea : public Gtk::GLArea
{
public:
    MyGLArea()
    {
        set_allowed_apis(Gdk::GLApi::GL);
        set_has_depth_buffer(true); // Explicitly request a depth buffer

        // Connect OpenGL event signals
        signal_realize().connect(sigc::mem_fun(*this, &MyGLArea::on_realize));
        signal_unrealize().connect(sigc::mem_fun(*this, &MyGLArea::on_unrealize), false);
        signal_render().connect(sigc::mem_fun(*this, &MyGLArea::on_render), false);
        signal_resize().connect(sigc::mem_fun(*this, &MyGLArea::on_resize)); // Connect resize signal

        // Initialize mouse control variables
        reset_camera_view_initial();
        m_highlighted_element = {NONE, -1}; // Initialize highlighted element

        // Connect mouse event signals
        auto motion_controller = Gtk::EventControllerMotion::create();
        motion_controller->signal_motion().connect(sigc::mem_fun(*this, &MyGLArea::on_mouse_motion));
        add_controller(motion_controller);

        auto click_controller_middle = Gtk::GestureClick::create();
        click_controller_middle->set_button(GDK_BUTTON_MIDDLE);
        click_controller_middle->signal_pressed().connect(
            [this](int n_press, double x, double y) { this->on_mouse_press(GDK_BUTTON_MIDDLE, x, y); }
        );
        click_controller_middle->signal_released().connect(
            [this](int n_press, double x, double y) { this->on_mouse_release(GDK_BUTTON_MIDDLE, x, y); }
        );
        add_controller(click_controller_middle);

        auto click_controller_right = Gtk::GestureClick::create();
        click_controller_right->set_button(GDK_BUTTON_SECONDARY);
        click_controller_right->signal_pressed().connect(
            [this](int n_press, double x, double y) { this->on_mouse_press(GDK_BUTTON_SECONDARY, x, y); }
        );
        click_controller_right->signal_released().connect(
            [this](int n_press, double x, double y) { this->on_mouse_release(GDK_BUTTON_SECONDARY, x, y); }
        );
        add_controller(click_controller_right);

        auto scroll_controller = Gtk::EventControllerScroll::create();
        scroll_controller->set_flags(Gtk::EventControllerScroll::Flags::VERTICAL);
        scroll_controller->signal_scroll().connect(
            [this](double dx, double dy) { return this->on_scroll(dx, dy); }, false
        );
        add_controller(scroll_controller);

        // Add key controller for Ctrl state
        auto key_controller = Gtk::EventControllerKey::create();
        key_controller->signal_key_pressed().connect(
            [this](guint keyval, guint keycode, Gdk::ModifierType state) -> bool { // Explicitly define return type as bool
                if (keyval == GDK_KEY_Control_L || keyval == GDK_KEY_Control_R) {
                    m_ctrl_pressed = true;
                }
                return false; // Don't stop propagation
            }, false // Explicitly pass 'false' for the 'after' parameter
        );
        key_controller->signal_key_released().connect(
            [this](guint keyval, guint keycode, Gdk::ModifierType state) { // Removed -> bool, returns void
                if (keyval == GDK_KEY_Control_L || keyval == GDK_KEY_Control_R) {
                    m_ctrl_pressed = false;
                }
                // No return statement needed for void
            }, false // Explicitly pass 'false' for the 'after' parameter
        );
        add_controller(key_controller);

        set_focusable(true);
        set_focus_on_click(true);
    }

    void load_step_file(const std::string& filepath)
    {
        try {
            make_current();
            throw_if_error();
            GL_CHECK_ERROR;

            std::cerr << "Loading STEP file: " << filepath << std::endl;
            // Corrected type: Dune3D::Result instead of Dune3D::STEPImporter::Result
            Dune3D::Result model_data = Dune3D::STEPImporter::import(filepath);

            if (model_data.faces.empty() && model_data.edges.empty()) {
                std::cerr << "Failed to load STEP file or no data found." << std::endl;
                // Optionally, load a default cube if STEP fails
                generate_cube_with_opencascade();
            } else {
                std::cerr << "STEP file loaded successfully. Faces: " << model_data.faces.size()
                          << ", Edges: " << model_data.edges.size() << std::endl;

                // Clear previous data
                m_vertex_data.clear();
                m_face_index_data.clear();
                m_edge_vertex_data.clear();
                m_edge_line_indices.clear();
                m_face_draw_info.clear(); // Clear draw info
                m_edge_draw_info.clear(); // Clear draw info
                m_face_surface_types.clear(); // Clear surface types

                // Populate face data
                GLuint current_base_index = 0;
                for (const auto& face : model_data.faces) {
                    size_t current_face_start_index = m_face_index_data.size(); // Store current size before adding indices
                    m_face_surface_types.push_back(face.surface_type); // Store the surface type for this logical face
                    for (const auto& vertex : face.vertices) {
                        m_vertex_data.push_back(vertex.x);
                        m_vertex_data.push_back(vertex.y);
                        m_vertex_data.push_back(vertex.z);
                        m_vertex_data.push_back(1.0f); // w component

                        m_vertex_data.push_back(face.color.r);
                        m_vertex_data.push_back(face.color.g);
                        m_vertex_data.push_back(face.color.b);
                        m_vertex_data.push_back(1.0f); // Alpha component
                    }
                    for (const auto& tri_indices : face.triangle_indices) {
                        m_face_index_data.push_back(current_base_index + tri_indices.x);
                        m_face_index_data.push_back(current_base_index + tri_indices.y);
                        m_face_index_data.push_back(current_base_index + tri_indices.z);
                    }
                    m_face_draw_info.emplace_back(current_face_start_index, m_face_index_data.size() - current_face_start_index);
                    current_base_index += static_cast<GLuint>(face.vertices.size());
                }

                // Populate edge data
                GLuint current_edge_base_index = 0;
                for (const auto& edge : model_data.edges) {
                    if (edge.size() < 2) continue; // Need at least two points for a line segment
                    size_t current_edge_start_index = m_edge_line_indices.size(); // Store current size before adding indices
                    for (const auto& vertex : edge) {
                        m_edge_vertex_data.push_back(vertex.x);
                        m_edge_vertex_data.push_back(vertex.y);
                        m_edge_vertex_data.push_back(vertex.z);
                        m_edge_vertex_data.push_back(1.0f); // w component
                    }
                    for (size_t i = 0; i < edge.size() - 1; ++i) {
                        m_edge_line_indices.push_back(current_edge_base_index + static_cast<GLuint>(i));
                        m_edge_line_indices.push_back(current_edge_base_index + static_cast<GLuint>(i + 1));
                    }
                    m_edge_draw_info.emplace_back(current_edge_start_index, m_edge_line_indices.size() - current_edge_start_index);
                    current_edge_base_index += static_cast<GLuint>(edge.size());
                }
            }

            calculate_bounding_box(); // Calculate bounding box after data is loaded
            init_buffers(); // Re-initialize buffers with new data
            reset_camera_view(); // Reset camera to fit the loaded model
            queue_draw();
        } catch (const Gdk::GLError& gle) {
            std::cerr << "GLArea error loading STEP: " << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error loading STEP file: " << e.what() << std::endl;
        }
    }

    // Public method to toggle culling
    void toggle_culling() {
        m_culling_enabled = !m_culling_enabled;
        std::cerr << "Culling " << (m_culling_enabled ? "enabled" : "disabled") << std::endl;
        queue_draw(); // Request redraw to apply the change
    }

    // Public method to toggle dynamic clipping plane
    void toggle_dynamic_clipping() {
        m_dynamic_clipping_enabled = !m_dynamic_clipping_enabled;
        std::cerr << "Dynamic Clipping Plane " << (m_dynamic_clipping_enabled ? "enabled" : "disabled") << std::endl;
        queue_draw(); // Request redraw to apply the change
    }

    // Public method to toggle projection mode
    void toggle_projection_mode() {
        m_is_perspective_projection = !m_is_perspective_projection;
        std::cerr << "Projection mode changed to " << (m_is_perspective_projection ? "Perspective" : "Orthographic") << std::endl;
        reset_camera_view(); // Reset camera view to adapt to new projection
        queue_draw(); // Request redraw to apply the change
    }

    // Public method to reset camera view to fit the model
    void reset_camera_view()
    {
        if (!m_has_model_data) {
            // Fallback to initial fixed view if no model data is loaded
            reset_camera_view_initial();
            return;
        }

        // Calculate center of the bounding box
        m_camera_target = m_model_center;
        m_rotation_pivot = m_model_center; // Set pivot to model center

        // Determine the largest dimension (X or Y) to fit the view, considering aspect ratio
        float aspect_ratio = (float)get_width() / (float)get_height();
        if (aspect_ratio == 0) aspect_ratio = 1.0f; // Prevent division by zero

        float target_view_height_from_y = m_model_dimensions.y;
        float target_view_height_from_x = m_model_dimensions.x / aspect_ratio;

        // Calculate base size for fitting the model, with padding
        float base_fit_size = std::max(target_view_height_from_y, target_view_height_from_x) * 1.25f;
        base_fit_size = std::max(0.1f, base_fit_size);

        // Update zoom levels based on the calculated base_fit_size
        m_ortho_view_size = base_fit_size;
        m_ortho_camera_distance = base_fit_size * 2.0f; // Maintain proportionality for orthographic camera distance

        const float FOV_DEG = 45.0f;
        m_perspective_camera_distance = (base_fit_size / 2.0f) / tan(glm::radians(FOV_DEG / 2.0f));

        // Preserve current camera viewing direction relative to target
        glm::vec3 current_camera_to_target_vec = m_camera_target - m_camera_position;
        glm::vec3 current_view_direction = glm::normalize(current_camera_to_target_vec);

        // If camera and target are at the same position, use a default direction
        if (glm::length(current_camera_to_target_vec) < 0.001f) {
            current_view_direction = glm::vec3(0.0f, 0.0f, -1.0f); // Default to looking down negative Z
        }

        // Update camera position based on the current projection mode
        float new_distance = m_is_perspective_projection ? m_perspective_camera_distance : m_ortho_camera_distance;
        m_camera_position = m_camera_target - current_view_direction * new_distance;

        // This part was commented out to prevent rotation reset on projection toggle.
        // If you want to reset rotation on general "reset camera view", uncomment this.
        // glm::quat rotX = glm::angleAxis(glm::radians(327.479f), glm::vec3(1.0f, 0.0f, 0.0f));
        // glm::quat rotY = glm::angleAxis(glm::radians(142.582f), glm::vec3(0.0f, 1.0f, 0.0f));
        // m_model_rotation = rotY * rotX; // Apply Y then X for consistent behavior

        queue_draw(); // Request redraw to apply the change
    }

    // Public method to set the view to "Front"
    void set_front_view() {
        if (!m_has_model_data) {
            reset_camera_view_initial();
            return;
        }

        m_model_rotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f); // Identity quaternion
        m_camera_target = m_model_center;
        m_rotation_pivot = m_model_center;

        float aspect_ratio = (float)get_width() / (float)get_height();
        if (aspect_ratio == 0) aspect_ratio = 1.0f;

        float target_view_height_from_y = m_model_dimensions.y;
        float target_view_height_from_x = m_model_dimensions.x / aspect_ratio;
        float base_fit_size = std::max(target_view_height_from_y, target_view_height_from_x) * 1.25f;
        base_fit_size = std::max(0.1f, base_fit_size);

        m_ortho_view_size = base_fit_size;
        m_perspective_camera_distance = (base_fit_size / 2.0f) / tan(glm::radians(45.0f / 2.0f));
        m_ortho_camera_distance = base_fit_size * 2.0f;

        float new_distance = m_is_perspective_projection ? m_perspective_camera_distance : m_ortho_camera_distance;
        m_camera_position = m_camera_target + glm::vec3(0.0f, 0.0f, new_distance);

        queue_draw();
    }

    // Public method to set the view to "Top"
    void set_top_view() {
        if (!m_has_model_data) {
            reset_camera_view_initial();
            return;
        }

        // Rotate model so that its top (positive Y) faces the camera (along world -Z)
        // This means rotating the model +90 degrees around the X-axis from its default orientation
        m_model_rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        m_camera_target = m_model_center;
        m_rotation_pivot = m_model_center;

        float aspect_ratio = (float)get_width() / (float)get_height();
        if (aspect_ratio == 0) aspect_ratio = 1.0f;

        // For top view, the relevant dimensions for fitting are X and Z
        float target_view_height_from_z = m_model_dimensions.z;
        float target_view_height_from_x = m_model_dimensions.x / aspect_ratio;
        float base_fit_size = std::max(target_view_height_from_z, target_view_height_from_x) * 1.25f;
        base_fit_size = std::max(0.1f, base_fit_size);

        m_ortho_view_size = base_fit_size;
        m_perspective_camera_distance = (base_fit_size / 2.0f) / tan(glm::radians(45.0f / 2.0f));
        m_ortho_camera_distance = base_fit_size * 2.0f;

        // Position the camera on the world's positive Z axis, looking towards negative Z
        float new_distance = m_is_perspective_projection ? m_perspective_camera_distance : m_ortho_camera_distance;
        m_camera_position = m_camera_target + glm::vec3(0.0f, 0.0f, new_distance);

        queue_draw();
    }

    // Public method to set the view to "Right"
    void set_right_view() {
        if (!m_has_model_data) {
            reset_camera_view_initial();
            return;
        }

        // Rotate model so that its right side (positive X) faces the camera (along world -Z)
        // This means rotating the model -90 degrees around the Y-axis from its default orientation
        m_model_rotation = glm::angleAxis(glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        m_camera_target = m_model_center;
        m_rotation_pivot = m_model_center;

        float aspect_ratio = (float)get_width() / (float)get_height();
        if (aspect_ratio == 0) aspect_ratio = 1.0f;

        // For right view, the relevant dimensions for fitting are Y and Z
        float target_view_height_from_y = m_model_dimensions.y;
        float target_view_height_from_z = m_model_dimensions.z / aspect_ratio;
        float base_fit_size = std::max(target_view_height_from_y, target_view_height_from_z) * 1.25f;
        base_fit_size = std::max(0.1f, base_fit_size);

        m_ortho_view_size = base_fit_size;
        m_perspective_camera_distance = (base_fit_size / 2.0f) / tan(glm::radians(45.0f / 2.0f));
        m_ortho_camera_distance = base_fit_size * 2.0f;

        // Position the camera on the world's positive Z axis, looking towards negative Z
        float new_distance = m_is_perspective_projection ? m_perspective_camera_distance : m_ortho_camera_distance;
        m_camera_position = m_camera_target + glm::vec3(0.0f, 0.0f, new_distance);

        queue_draw();
    }

    // Public method to toggle normal correction tint
    void toggle_normal_correction_tint() {
        m_apply_normal_correction_tint = !m_apply_normal_correction_tint;
        std::cerr << "Normal correction tint " << (m_apply_normal_correction_tint ? "enabled" : "disabled") << std::endl;
        queue_draw();
    }


protected:
    void on_realize() override
    {
        Gtk::GLArea::on_realize();
        try {
            make_current();
            throw_if_error();
            std::cerr << "GLArea realized and made current." << std::endl;
            GL_CHECK_ERROR;

            generate_cube_with_opencascade(); // Initial cube for demonstration
            std::cerr << "Initial cube data generated with OpenCASCADE." << std::endl;
            calculate_bounding_box(); // Calculate bounding box for the initial cube
            init_buffers();        // Initialize OpenGL buffers
            std::cerr << "Buffers initialized." << std::endl;
            init_shaders();        // Compile and link shaders
            std::cerr << "Shaders initialized." << std::endl;

            glEnable(GL_DEPTH_TEST); // Enable depth testing
            // Initial culling state based on m_culling_enabled
            if (m_culling_enabled) {
                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);
            } else {
                glDisable(GL_CULL_FACE);
            }

            glClearColor(0.6f, 0.7f, 0.9f, 1.0f); // Light blue background
            std::cerr << "Global OpenGL states set." << std::endl;
            GL_CHECK_ERROR;

        } catch (const Gdk::GLError& gle) {
            std::cerr << "GLArea error on realize: " << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
        }
    }

    void on_unrealize() override
    {
        try {
            make_current();
            throw_if_error();
            std::cerr << "GLArea unrealizing." << std::endl;
            GL_CHECK_ERROR;

            // Delete OpenGL buffers and program to free resources
            glDeleteBuffers(1, &m_position_buffer);
            glDeleteBuffers(1, &m_index_buffer);
            glDeleteBuffers(1, &m_edge_position_buffer);
            glDeleteBuffers(1, &m_edge_index_buffer);
            glDeleteProgram(m_program);

            std::cerr << "OpenGL resources deleted." << std::endl;
            GL_CHECK_ERROR;

        } catch (const Gdk::GLError& gle) {
            std::cerr << "GLArea error on unrealize: " << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
        }

        Gtk::GLArea::on_unrealize();
    }

    // Method to handle resize
    void on_resize(int width, int height) override
    {
        std::cerr << "on_resize called with dimensions: " << width << "x" << height << std::endl;
        try {
            make_current();
            throw_if_error();
            GL_CHECK_ERROR;

            // Set the viewport
            glViewport(0, 0, width, height);
            GL_CHECK_ERROR;

            // Recalculate camera position to re-fit the model to the new window size
            reset_camera_view();

        } catch (const Gdk::GLError& gle) {
            std::cerr << "GLArea error on resize: " << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
        }
    }

    // Converts 2D mouse coordinates to a 3D point on a virtual trackball sphere
    glm::vec3 map_to_trackball(float x, float y, float width, float height)
    {
        glm::vec3 v;
        // Normalize coordinates to range [-1, 1]
        v.x = (2.0f * x - width) / width;
        v.y = (height - 2.0f * y) / height; // Y is inverted for OpenGL
        v.z = 0.0f;

        float op_sq = v.x * v.x + v.y * v.y;
        if (op_sq <= 1.0f) {
            v.z = std::sqrt(1.0f - op_sq);
        } else {
            // If outside the circle, normalize and project onto the sphere's edge
            // This case should ideally not happen if the trackball is well-defined
            v = glm::normalize(v);
        }
        return v;
    }

    // Mouse motion handler
    void on_mouse_motion(double x, double double_y)
    {
        if (m_is_rotating) {
            glm::vec3 current_pos_3d = map_to_trackball(static_cast<float>(x), static_cast<float>(double_y),
                                                        static_cast<float>(get_width()), static_cast<float>(get_height()));
            glm::vec3 last_pos_3d = map_to_trackball(m_last_mouse_pos.x, m_last_mouse_pos.y,
                                                     static_cast<float>(get_width()), static_cast<float>(get_height()));

            // Calculate rotation quaternion
            glm::quat new_rotation = glm::rotation(last_pos_3d, current_pos_3d);

            // Apply new rotation to the current model rotation
            m_model_rotation = new_rotation * m_model_rotation;
            m_model_rotation = glm::normalize(m_model_rotation); // Normalize for stability

            queue_draw();
        } else if (m_is_panning) {
            float delta_x = static_cast<float>(x - m_last_mouse_pos.x);
            float delta_y = static_cast<float>(double_y - m_last_mouse_pos.y);

            // Pan speed should be relative to the current orthographic view size
            float current_pan_scale = m_ortho_view_size / (float)get_height(); // Scale by view height

            glm::vec3 camera_direction = glm::normalize(m_camera_target - m_camera_position);
            glm::vec3 world_up = glm::vec3(0.0f, 1.0f, 0.0f);
            glm::vec3 camera_right = glm::normalize(glm::cross(camera_direction, world_up));
            glm::vec3 camera_up = glm::normalize(glm::cross(camera_right, camera_direction));

            m_camera_position -= camera_right * delta_x * m_pan_speed * current_pan_scale;
            m_camera_target -= camera_right * delta_x * m_pan_speed * current_pan_scale;

            m_camera_position += camera_up * delta_y * m_pan_speed * current_pan_scale;
            m_camera_target += camera_up * delta_y * m_pan_speed * current_pan_scale;

            queue_draw();
        }
        m_last_mouse_pos = glm::vec2(x, double_y);

        // Picking logic when not rotating or panning
        if (!m_is_rotating && !m_is_panning) {
            // Get current view and projection matrices
            glm::mat4 projection = compute_projection_matrix();
            glm::mat4 view = glm::lookAt(m_camera_position, m_camera_target, glm::vec3(0.0f, 1.0f, 0.0f));

            // Define the full model transformation matrix for consistency with rendering
            // This should match the model matrix calculation in compute_mvp()
            glm::mat4 model_transform = glm::translate(glm::mat4(1.0f), m_rotation_pivot) *
                                        glm::mat4_cast(m_model_rotation) *
                                        glm::translate(glm::mat4(1.0f), -m_rotation_pivot);

            // Corrected: Ray should be calculated in world space (after view and projection, but BEFORE model_transform)
            // The vertices for picking are already transformed by model_transform, so the ray needs to be in the same space.
            glm::mat4 inv_view_proj = glm::inverse(projection * view);

            // Ray from mouse position (in NDC)
            float ndc_x = (2.0f * static_cast<float>(x)) / get_width() - 1.0f;
            float ndc_y = 1.0f - (2.0f * static_cast<float>(double_y)) / get_height();

            // Create ray in world space
            glm::vec4 ray_clip_near = glm::vec4(ndc_x, ndc_y, -1.0f, 1.0f);
            glm::vec4 ray_clip_far = glm::vec4(ndc_x, ndc_y, 1.0f, 1.0f);

            glm::vec4 ray_world_near = inv_view_proj * ray_clip_near;
            glm::vec4 ray_world_far = inv_view_proj * ray_clip_far;

            glm::vec3 ray_origin = glm::vec3(ray_world_near) / ray_world_near.w;
            glm::vec3 ray_direction = glm::normalize(glm::vec3(ray_world_far) / ray_world_far.w - ray_origin);

            HighlightedElement new_highlight = {NONE, -1};
            float min_vertex_dist_sq = std::numeric_limits<float>::max();
            int picked_vertex_id = -1;

            float min_edge_dist_sq = std::numeric_limits<float>::max();
            int picked_edge_id = -1;

            float min_face_t = std::numeric_limits<float>::max();
            int picked_face_id = -1;

            // 1. Pick Faces (highest priority)
            for (size_t i = 0; i < m_face_draw_info.size(); ++i) { // Iterate through each logical face
                const auto& draw_info = m_face_draw_info[i];
                for (size_t j = 0; j < draw_info.second; j += 3) {
                    GLuint idx0 = m_face_index_data[draw_info.first + j];
                    GLuint idx1 = m_face_index_data[draw_info.first + j + 1];
                    GLuint idx2 = m_face_index_data[draw_info.first + j + 2];

                    glm::vec3 v0(m_vertex_data[idx0 * 8], m_vertex_data[idx0 * 8 + 1], m_vertex_data[idx0 * 8 + 2]);
                    glm::vec3 v1(m_vertex_data[idx1 * 8], m_vertex_data[idx1 * 8 + 1], m_vertex_data[idx1 * 8 + 2]);
                    glm::vec3 v2(m_vertex_data[idx2 * 8], m_vertex_data[idx2 * 8 + 1], m_vertex_data[idx2 * 8 + 2]);

                    glm::vec3 transformed_v0 = glm::vec3(model_transform * glm::vec4(v0, 1.0f));
                    glm::vec3 transformed_v1 = glm::vec3(model_transform * glm::vec4(v1, 1.0f));
                    glm::vec3 transformed_v2 = glm::vec3(model_transform * glm::vec4(v2, 1.0f));

                    float t;
                    if (ray_triangle_intersect(ray_origin, ray_direction, transformed_v0, transformed_v1, transformed_v2, t)) {
                        if (t > 0 && t < min_face_t) { // Ensure intersection is in front of the ray origin
                            min_face_t = t;
                            picked_face_id = static_cast<int>(i);
                        }
                    }
                }
            }

            // 2. Pick Edges (medium priority)
            // Increased threshold for easier picking of edges
            float current_edge_picking_threshold_sq = 1.5f * 1.5f; // Increased from 0.8*0.8
            for (size_t i = 0; i < m_edge_draw_info.size(); ++i) { // Iterate through each logical edge
                const auto& draw_info = m_edge_draw_info[i];
                for (size_t j = 0; j < draw_info.second; j += 2) {
                    GLuint idx0 = m_edge_line_indices[draw_info.first + j];
                    GLuint idx1 = m_edge_line_indices[draw_info.first + j + 1];

                    glm::vec3 p0(m_edge_vertex_data[idx0 * 4], m_edge_vertex_data[idx0 * 4 + 1], m_edge_vertex_data[idx0 * 4 + 2]);
                    glm::vec3 p1(m_edge_vertex_data[idx1 * 4], m_edge_vertex_data[idx1 * 4 + 1], m_edge_vertex_data[idx1 * 4 + 2]);

                    glm::vec3 transformed_p0 = glm::vec3(model_transform * glm::vec4(p0, 1.0f));
                    glm::vec3 transformed_p1 = glm::vec3(model_transform * glm::vec4(p1, 1.0f));

                    glm::vec3 segment_dir = transformed_p1 - transformed_p0;
                    float segment_length_sq = glm::dot(segment_dir, segment_dir);

                    float t_ray = 0.0f;
                    float t_segment = 0.0f;

                    if (segment_length_sq > 0.0f) {
                        glm::vec3 w0 = ray_origin - transformed_p0;
                        float a = glm::dot(ray_direction, ray_direction);
                        float b = glm::dot(ray_direction, segment_dir);
                        float c = glm::dot(segment_dir, segment_dir);
                        float d = glm::dot(ray_direction, w0);
                        float e = glm::dot(segment_dir, w0);

                        float denom = a * c - b * b;
                        if (denom < 0.00001f) { // Parallel or nearly parallel
                            t_segment = glm::dot(ray_origin - transformed_p0, segment_dir) / segment_length_sq;
                            t_segment = glm::clamp(t_segment, 0.0f, 1.0f);
                            t_ray = 0.0f; // Arbitrary, since parallel
                        } else {
                            t_ray = (b * e - c * d) / denom;
                            t_segment = (a * e - b * d) / denom;
                            t_segment = glm::clamp(t_segment, 0.0f, 1.0f);
                        }
                    }

                    glm::vec3 closest_point_on_ray = ray_origin + t_ray * ray_direction;
                    glm::vec3 closest_point_on_segment = transformed_p0 + t_segment * segment_dir;

                    float dist_sq = glm::distance2(closest_point_on_ray, closest_point_on_segment);

                    if (dist_sq < current_edge_picking_threshold_sq && dist_sq < min_edge_dist_sq) {
                        min_edge_dist_sq = dist_sq;
                        picked_edge_id = static_cast<int>(i);
                    }
                }
            }

            // 3. Pick Vertices (lowest priority)
            // Increased threshold for easier picking of vertices
            float current_vertex_picking_threshold_sq = 1.0f * 1.0f; // Increased from 0.5*0.5
            for (size_t i = 0; i < m_vertex_data.size() / 8; ++i) { // Divide by 8 because each vertex has 8 floats
                glm::vec3 v_pos(m_vertex_data[i * 8], m_vertex_data[i * 8 + 1], m_vertex_data[i * 8 + 2]);
                // transformed_v is the vertex in its final world position after user interaction
                glm::vec3 transformed_v = glm::vec3(model_transform * glm::vec4(v_pos, 1.0f));

                glm::vec3 vec_to_vertex = transformed_v - ray_origin;
                float t_on_ray = glm::dot(vec_to_vertex, ray_direction);
                glm::vec3 closest_point_on_ray_to_vertex = ray_origin + t_on_ray * ray_direction;

                float dist_sq = glm::distance2(transformed_v, closest_point_on_ray_to_vertex);

                if (dist_sq < current_vertex_picking_threshold_sq && dist_sq < min_vertex_dist_sq) {
                    min_vertex_dist_sq = dist_sq;
                    picked_vertex_id = static_cast<int>(i);
                }
            }


            // Determine final highlight based on hierarchy: Face > Edge > Vertex
            if (picked_face_id != -1) {
                new_highlight = {FACE, picked_face_id};
            } else if (picked_edge_id != -1) {
                new_highlight = {EDGE, picked_edge_id};
            } else if (picked_vertex_id != -1) {
                new_highlight = {VERTEX, picked_vertex_id};
                // Diagnostic output: Mouse coordinates (now only when a vertex is highlighted)
                std::cerr << "Mouse (Screen): (" << x << ", " << double_y << ")" << std::endl;
                // Diagnostic output: Highlighted vertex coordinates
                glm::vec3 v_pos(m_vertex_data[picked_vertex_id * 8], m_vertex_data[picked_vertex_id * 8 + 1], m_vertex_data[picked_vertex_id * 8 + 2]);
                glm::vec3 transformed_v = glm::vec3(model_transform * glm::vec4(v_pos, 1.0f));
                std::cerr << "Highlighted Vertex (World): (" << transformed_v.x << ", " << transformed_v.y << ", " << transformed_v.z << ")" << std::endl;

                // Convert highlighted vertex world coordinates to screen coordinates
                // Use the already transformed world vertex (transformed_v)
                // and apply only view and projection matrices.
                glm::vec4 clip_pos = projection * view * glm::vec4(transformed_v, 1.0f);
                glm::vec3 ndc_pos = glm::vec3(clip_pos) / clip_pos.w;

                int screen_width = get_width();
                int screen_height = get_height();
                // Convert NDC to screen coordinates, matching the mouse coordinate system
                float highlighted_screen_x = (ndc_pos.x + 1.0f) * screen_width / 2.0f;
                float highlighted_screen_y = (1.0f - ndc_pos.y) * screen_height / 2.0f;

                std::cerr << "Highlighted Vertex (Screen): (" << highlighted_screen_x << ", " << highlighted_screen_y << ")" << std::endl;
            }


            if (new_highlight.type != m_highlighted_element.type || new_highlight.id != m_highlighted_element.id) {
                m_highlighted_element = new_highlight;
                queue_draw();
            }
        }
    }

    // Mouse button press handler
    void on_mouse_press(int button, double x, double y)
    {
        m_last_mouse_pos = glm::vec2(x, y);

        if (button == GDK_BUTTON_MIDDLE) { // Middle button for rotation
            m_middle_button_pressed = true;
        } else if (button == GDK_BUTTON_SECONDARY) { // Right button for panning / setting pivot
            m_right_button_pressed = true;
        }

        // Clear highlight when any interaction starts
        m_highlighted_element = {NONE, -1};
        queue_draw(); // Request redraw to clear highlight

        // Logic for setting/resetting rotation pivot
        if (m_middle_button_pressed && m_right_button_pressed) {
            if (m_ctrl_pressed) {
                // Reset pivot to model center
                m_rotation_pivot = m_model_center;
                std::cerr << "Rotation pivot reset to model center: (" << m_rotation_pivot.x << ", " << m_rotation_pivot.y << ", " << m_rotation_pivot.z << ")" << std::endl;
            } else {
                // Set new pivot using raycasting
                set_rotation_pivot_from_mouse(x, y);
            }
            // After setting pivot, still assume rotation mode
            m_is_rotating = true;
            m_is_panning = false;
        } else if (m_middle_button_pressed) {
            m_is_rotating = true;
            m_is_panning = false;
        } else if (m_right_button_pressed) {
            m_is_panning = true;
            m_is_rotating = false;
        }
        grab_focus(); // Get focus to receive events
        queue_draw();
    }

    // Mouse button release handler
    void on_mouse_release(int button, double x, double y)
    {
        if (button == GDK_BUTTON_MIDDLE) {
            m_middle_button_pressed = false;
        } else if (button == GDK_BUTTON_SECONDARY) {
            m_right_button_pressed = false;
        }

        // If both buttons are released, stop all actions
        if (!m_middle_button_pressed && !m_right_button_pressed) {
            m_is_rotating = false;
            m_is_panning = false;
            // Re-evaluate highlight after interaction ends
            on_mouse_motion(x, y); // Trigger picking for the current mouse position
        }
        // If one button is still pressed, continue its action
        else if (m_middle_button_pressed) {
            m_is_rotating = true;
            m_is_panning = false;
        } else if (m_right_button_pressed) {
            m_is_panning = true;
            m_is_rotating = false;
        }
        queue_draw();
    }

    // Converts 2D screen coordinates to a 3D point in world space on the depth plane
    glm::vec3 unproject_mouse_to_world(double mouse_x, double mouse_y, float depth_z) {
        glm::mat4 projection = compute_projection_matrix();
        glm::mat4 view = glm::lookAt(m_camera_position, m_camera_target, glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 inv_mvp = glm::inverse(projection * view); // This is correct for unprojecting to world space

        // Normalized Device Coordinates (NDC)
        // Corrected: No inversion for ndc_x here.
        float ndc_x = (2.0f * static_cast<float>(mouse_x)) / get_width() - 1.0f;
        float ndc_y = 1.0f - (2.0f * static_cast<float>(mouse_y)) / get_height(); // Y is inverted for OpenGL
        float ndc_z = depth_z; // Depth, usually in range [0, 1] for z_near and z_far

        glm::vec4 screen_pos = glm::vec4(ndc_x, ndc_y, ndc_z, 1.0f);
        glm::vec4 world_pos = inv_mvp * screen_pos;
        return glm::vec3(world_pos) / world_pos.w;
    }

    // Mouse scroll wheel handler (zoom)
    gboolean on_scroll(double dx, double double_y)
    {
        // Clear highlight when zooming starts
        m_highlighted_element = {NONE, -1};

        // Get current mouse position (use last known position from on_mouse_motion)
        double mouse_x = m_last_mouse_pos.x;
        double mouse_y = m_last_mouse_pos.y;

        float zoom_factor = 1.0f;
        if (double_y > 0) { // Scroll up (zoom in)
            zoom_factor = 1.0f / 1.1f; // Zoom in by 10%
        } else if (double_y < 0) { // Scroll down (zoom out)
            zoom_factor = 1.1f; // Zoom out by 10%
        }

        // Calculate the point in world space that is currently under the mouse cursor
        // This point will be used as the center of the zoom operation.
        glm::vec3 zoom_center_world = unproject_mouse_to_world(mouse_x, mouse_y, 0.5f);

        if (m_is_perspective_projection) {
            // Adjust perspective camera distance directly
            m_perspective_camera_distance *= zoom_factor;
            // Clamp to a much smaller minimum to allow very close zoom
            m_perspective_camera_distance = glm::clamp(m_perspective_camera_distance, 0.001f, 10000.0f);

            // Recalculate camera position based on new distance, preserving view direction
            glm::vec3 current_camera_to_target_vec = m_camera_target - m_camera_position;
            glm::vec3 view_direction = glm::normalize(current_camera_to_target_vec);
            m_camera_position = m_camera_target - view_direction * m_perspective_camera_distance;

        } else {
            // Adjust orthographic view size
            m_ortho_view_size *= zoom_factor;
            m_ortho_view_size = glm::clamp(m_ortho_view_size, 0.001f, 10000.0f); // Also clamp to smaller minimum

            // Also adjust orthographic camera distance to simulate camera movement
            m_ortho_camera_distance *= zoom_factor;
            m_ortho_camera_distance = glm::clamp(m_ortho_camera_distance, 0.001f, 10000.0f); // Also clamp to smaller minimum

            // Recalculate camera position based on new distance, preserving view direction
            glm::vec3 current_camera_to_target_vec = m_camera_target - m_camera_position;
            glm::vec3 view_direction = glm::normalize(current_camera_to_target_vec);
            m_camera_position = m_camera_target - view_direction * m_ortho_camera_distance;
        }

        // After adjusting zoom level/distance, the projection matrix changes.
        // We need to re-unproject the mouse position to find its new world coordinate.
        glm::vec3 new_zoom_center_world = unproject_mouse_to_world(mouse_x, mouse_y, 0.5f);

        // Calculate the displacement needed to keep the original zoom_center_world fixed under the mouse.
        // This effectively pans the camera to compensate for the zoom, including perpendicular movement.
        glm::vec3 pan_correction = zoom_center_world - new_zoom_center_world;

        // Apply the displacement to both camera position and target
        m_camera_position += pan_correction;
        m_camera_target += pan_correction;

        queue_draw();
        return TRUE;
    }


    bool on_render(const Glib::RefPtr<Gdk::GLContext>& context) override
    {
        try {
            make_current();
            throw_if_error();
            GL_CHECK_ERROR;

            // Get viewport dimensions (for MVP matrix)
            int width = get_width();
            int height = get_height();

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            GL_CHECK_ERROR;
            draw_object(); // Render
            GL_CHECK_ERROR;
            glFlush(); // Ensure rendering is complete
            GL_CHECK_ERROR;
            return true;
        } catch (const Gdk::GLError& gle) {
            std::cerr << "GLArea error on render: " << gle.domain() << "-" << gle.code() << "-" << gle.what() << std::endl;
            return false;
        }
    }

private:
    GLuint m_position_buffer = 0;    // VBO for face positions and colors
    GLuint m_index_buffer = 0;       // EBO for face indices
    GLuint m_edge_position_buffer = 0; // VBO for edge positions
    GLuint m_edge_index_buffer = 0;    // EBO for edge indices
    GLuint m_program = 0;            // Shader program ID
    GLint m_mvp_location = -1;       // Location of MVP matrix uniform
    GLint m_is_line_pass_location = -1; // Location of uniform to determine pass (faces/edges)
    GLint m_line_color_uniform_location = -1; // Location of uniform for line color

    // Variables for mouse and camera control
    glm::quat m_model_rotation; // Model rotation using quaternion
    glm::vec3 m_rotation_pivot; // New rotation point
    glm::vec2 m_last_mouse_pos;
    bool m_is_rotating = false;
    bool m_is_panning = false;
    bool m_middle_button_pressed = false; // Track middle button state
    bool m_right_button_pressed = false;  // Track right button state
    bool m_ctrl_pressed = false;          // Track Ctrl key state

    glm::vec3 m_camera_position;
    glm::vec3 m_camera_target;
    float m_rotation_speed = 0.7f; // This speed is now less direct with trackball
    float m_pan_speed = 1.0f;
    float m_zoom_speed = 0.4f;
    bool m_culling_enabled = true; // New member to control culling state
    bool m_dynamic_clipping_enabled = false; // New member for dynamic clipping plane
    bool m_is_perspective_projection = false; // New member for projection switching

    // Separate zoom variables for orthographic and perspective projection
    float m_ortho_view_size = 2.0f; // Initial orthographic view size (height of frustum)
    float m_ortho_camera_distance = 2.0f; // New: Distance of camera from target in ortho mode
    float m_perspective_camera_distance = 2.0f; // Initial camera distance from target for perspective

    std::vector<GLfloat> m_vertex_data;      // Vertex data (position + color) for faces
    std::vector<GLuint> m_face_index_data;  // Indices for drawing faces (triangles)
    std::vector<GLfloat> m_edge_vertex_data; // Vertex data (position only) for edges
    std::vector<GLuint> m_edge_line_indices;  // Indices for drawing edges (lines)

    // Data for efficient drawing of highlighted elements
    std::vector<std::pair<size_t, size_t>> m_face_draw_info; // {start_index_in_m_face_index_data, count}
    std::vector<std::pair<size_t, size_t>> m_edge_draw_info; // {start_index_in_m_edge_line_indices, count}

    // Bounding box variables
    glm::vec3 m_min_bounds = glm::vec3(std::numeric_limits<float>::max());
    glm::vec3 m_max_bounds = glm::vec3(std::numeric_limits<float>::lowest());
    glm::vec3 m_model_center = glm::vec3(0.0f);
    glm::vec3 m_model_dimensions = glm::vec3(0.0f);
    bool m_has_model_data = false; // To check if bounding box is valid

    // Variables for highlighting
    HighlightedElement m_highlighted_element;
    std::vector<GeomAbs_SurfaceType> m_face_surface_types; // New: Store surface type for each logical face

    // New members for normal correction tint
    bool m_apply_normal_correction_tint = false;
    glm::vec4 m_normal_correction_tint_color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f); // Darker blue-grey tint
    GLint m_apply_normal_correction_tint_location = -1;
    GLint m_normal_correction_tint_color_location = -1;


    // Declare generate_cube_with_opencascade here
    void generate_cube_with_opencascade();

    // Resets camera view to a fixed initial position (used before model is loaded)
    void reset_camera_view_initial()
    {
        m_camera_target = glm::vec3(0.0f, 0.0f, 0.0f);
        m_rotation_pivot = glm::vec3(0.0f, 0.0f, 0.0f); // Initial cube is at origin
        // Initial rotation for quaternion (corresponds to old Euler angles)
        glm::quat rotX = glm::angleAxis(glm::radians(327.479f), glm::vec3(1.0f, 0.0f, 0.0f));
        glm::quat rotY = glm::angleAxis(glm::radians(142.582f), glm::vec3(0.0f, 1.0f, 0.0f));
        m_model_rotation = rotY * rotX;

        // Initial values for zoom levels based on a default cube size
        const float initial_cube_size = 50.0f;
        m_ortho_view_size = initial_cube_size * 1.25f; // For orthographic, height of the frustum
        m_ortho_camera_distance = initial_cube_size * 2.0f; // A reasonable initial distance for ortho camera
        const float FOV_DEG = 45.0f;
        // For perspective, calculate initial distance to fit the cube within the FOV
        m_perspective_camera_distance = (initial_cube_size * 1.25f / 2.0f) / tan(glm::radians(FOV_DEG / 2.0f));

        m_has_model_data = false; // No model data loaded yet for initial state

        // Set initial camera position along the world Z axis (looking at origin)
        // Use the appropriate distance based on the initial projection mode
        float initial_distance = m_is_perspective_projection ? m_perspective_camera_distance : m_ortho_camera_distance;
        m_camera_position = m_camera_target + glm::vec3(0.0f, 0.0f, 1.0f) * initial_distance;
    }

    // Calculates the bounding box of the loaded model
    void calculate_bounding_box()
    {
        m_min_bounds = glm::vec3(std::numeric_limits<float>::max());
        m_max_bounds = glm::vec3(std::numeric_limits<float>::lowest());
        m_has_model_data = false;

        // Iterate through face vertices
        for (size_t i = 0; i < m_vertex_data.size(); i += 8) { // 8 floats per vertex (pos + color)
            glm::vec3 pos(m_vertex_data[i], m_vertex_data[i+1], m_vertex_data[i+2]);
            m_min_bounds = glm::min(m_min_bounds, pos);
            m_max_bounds = glm::max(m_max_bounds, pos);
            m_has_model_data = true;
        }

        // Iterate through edge vertices
        for (size_t i = 0; i < m_edge_vertex_data.size(); i += 4) { // 4 floats per vertex (pos only)
            glm::vec3 pos(m_edge_vertex_data[i], m_edge_vertex_data[i+1], m_edge_vertex_data[i+2]);
            m_min_bounds = glm::min(m_min_bounds, pos);
            m_max_bounds = glm::max(m_max_bounds, pos);
            m_has_model_data = true;
        }

        if (m_has_model_data) {
            m_model_center = (m_min_bounds + m_max_bounds) / 2.0f;
            m_model_dimensions = m_max_bounds - m_min_bounds;
            std::cerr << "Bounding Box: Min(" << m_min_bounds.x << ", " << m_min_bounds.y << ", " << m_min_bounds.z
                      << "), Max(" << m_max_bounds.x << ", " << m_max_bounds.y << ", " << m_max_bounds.z
                      << "), Center(" << m_model_center.x << ", " << m_model_center.y << ", " << m_model_center.z
                      << "), Dimensions(" << m_model_dimensions.x << ", " << m_model_dimensions.y << ", " << m_model_dimensions.z << ")" << std::endl;
        } else {
            m_model_center = glm::vec3(0.0f);
            m_model_dimensions = glm::vec3(0.0f);
            std::cerr << "No model data to calculate bounding box." << std::endl;
        }
    }

    // Initializes OpenGL buffer objects (VBOs and EBOs)
    void init_buffers()
    {
        GLuint vao;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        // Delete existing buffers if they exist
        if (m_position_buffer != 0) glDeleteBuffers(1, &m_position_buffer);
        if (m_index_buffer != 0) glDeleteBuffers(1, &m_index_buffer);
        if (m_edge_position_buffer != 0) glDeleteBuffers(1, &m_edge_position_buffer);
        if (m_edge_index_buffer != 0) glDeleteBuffers(1, &m_edge_index_buffer);

        // Generate and bind buffers for faces
        glGenBuffers(1, &m_position_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, m_position_buffer);
        glBufferData(GL_ARRAY_BUFFER, m_vertex_data.size() * sizeof(GLfloat), m_vertex_data.data(), GL_STATIC_DRAW);

        glGenBuffers(1, &m_index_buffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_face_index_data.size() * sizeof(GLuint), m_face_index_data.data(), GL_STATIC_DRAW);

        // Generate and bind buffers for edges
        glGenBuffers(1, &m_edge_position_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, m_edge_position_buffer);
        glBufferData(GL_ARRAY_BUFFER, m_edge_vertex_data.size() * sizeof(GLfloat), m_edge_vertex_data.data(), GL_STATIC_DRAW);

        glGenBuffers(1, &m_edge_index_buffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_edge_index_buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_edge_line_indices.size() * sizeof(GLuint), m_edge_line_indices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    // Creates and compiles a shader
    GLuint create_shader(GLenum type, const char* src)
    {
        GLuint shader = glCreateShader(type);
        glShaderSource(shader, 1, &src, NULL);
        glCompileShader(shader);

        GLint status;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
        if (status == GL_FALSE) {
            GLint log_len;
            glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_len);
            std::vector<char> buffer(log_len + 1);
            glGetShaderInfoLog(shader, log_len, NULL, buffer.data());
            std::cerr << "Shader compilation failed (" << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << "):\n" << buffer.data() << std::endl;
            glDeleteShader(shader);
            return 0;
        }
        return shader;
    }

    // Initializes and links the shader program
    void init_shaders()
    {
        GLuint vertex_shader = create_shader(GL_VERTEX_SHADER, vertex_shader_src);
        GLuint fragment_shader = create_shader(GL_FRAGMENT_SHADER, fragment_shader_src);

        if (vertex_shader == 0 || fragment_shader == 0) {
            m_program = 0;
            return;
        }

        m_program = glCreateProgram();
        glAttachShader(m_program, vertex_shader);
        glAttachShader(m_program, fragment_shader);

        glBindAttribLocation(m_program, 0, "in_position");
        glBindAttribLocation(m_program, 1, "in_color");

        glLinkProgram(m_program);

        GLint status;
        glGetProgramiv(m_program, GL_LINK_STATUS, &status);
        if (status == GL_FALSE) {
            GLint log_len;
            glGetProgramiv(m_program, GL_INFO_LOG_LENGTH, &log_len);
            std::vector<char> buffer(log_len + 1);
            glGetProgramInfoLog(m_program, log_len, NULL, buffer.data());
            std::cerr << "Program linking failed:\n" << buffer.data() << std::endl;
            glDeleteProgram(m_program);
            m_program = 0;
            goto out;
        }

        m_mvp_location = glGetUniformLocation(m_program, "mvp");
        m_is_line_pass_location = glGetUniformLocation(m_program, "is_line_pass");
        m_line_color_uniform_location = glGetUniformLocation(m_program, "line_color_uniform");
        m_apply_normal_correction_tint_location = glGetUniformLocation(m_program, "apply_normal_correction_tint");
        m_normal_correction_tint_color_location = glGetUniformLocation(m_program, "normal_correction_tint_color");


    out:
        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
    }

    // Helper to get only the projection matrix
    glm::mat4 compute_projection_matrix() {
        float aspect_ratio = (float)get_width() / (float)get_height();
        if (aspect_ratio == 0) aspect_ratio = 1.0f;

        // Calculate the Z-range of the *transformed* model in camera space.
        glm::mat4 model_transform = glm::translate(glm::mat4(1.0f), m_rotation_pivot) *
                                    glm::mat4_cast(m_model_rotation) *
                                    glm::translate(glm::mat4(1.0f), -m_rotation_pivot);
        glm::mat4 view = glm::lookAt(m_camera_position, m_camera_target, glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 model_view = view * model_transform;

        float min_z_view_space = std::numeric_limits<float>::max(); // Most negative Z (farthest)
        float max_z_view_space = std::numeric_limits<float>::lowest(); // Least negative Z (closest)

        // Iterate through the bounding box corners (8 corners) in model space
        glm::vec3 corners[8];
        corners[0] = glm::vec3(m_min_bounds.x, m_min_bounds.y, m_min_bounds.z);
        corners[1] = glm::vec3(m_max_bounds.x, m_min_bounds.y, m_min_bounds.z);
        corners[2] = glm::vec3(m_min_bounds.x, m_max_bounds.y, m_min_bounds.z);
        corners[3] = glm::vec3(m_min_bounds.x, m_min_bounds.y, m_max_bounds.z);
        corners[4] = glm::vec3(m_max_bounds.x, m_max_bounds.y, m_min_bounds.z);
        corners[5] = glm::vec3(m_max_bounds.x, m_min_bounds.y, m_max_bounds.z);
        corners[6] = glm::vec3(m_min_bounds.x, m_max_bounds.y, m_max_bounds.z);
        corners[7] = glm::vec3(m_max_bounds.x, m_max_bounds.y, m_max_bounds.z);

        for (int i = 0; i < 8; ++i) {
            glm::vec4 transformed_corner = model_view * glm::vec4(corners[i], 1.0f);
            min_z_view_space = glm::min(min_z_view_space, transformed_corner.z);
            max_z_view_space = glm::max(max_z_view_space, transformed_corner.z);
        }

        float z_near_dist, z_far_dist;

        // Calculate the actual distances of the model's bounding box from the camera
        // (these will be positive values)
        float model_closest_dist = -max_z_view_space;
        float model_furthest_dist = -min_z_view_space;

        // Ensure minimum positive distance for near plane
        model_closest_dist = glm::max(0.001f, model_closest_dist); // Changed from 0.1f to 0.001f


        if (m_is_perspective_projection) {
            const float FOV_RAD = glm::radians(45.0f); // Fixed comfortable FOV for perspective

            if (m_dynamic_clipping_enabled) {
                const float FIXED_CLIP_DISTANCE = 60.0f; // 60mm from camera
                z_near_dist = FIXED_CLIP_DISTANCE;
                // z_far_dist should extend beyond the furthest part of the model that is visible
                // from the fixed near plane.
                z_far_dist = glm::max(z_near_dist + 0.1f, model_furthest_dist + 1.0f); // Add a small buffer
            } else {
                // Non-dynamic perspective clipping: encompass the whole model
                z_near_dist = glm::max(0.001f, model_closest_dist - 1.0f); // Extend slightly in front, changed from 0.1f
                z_far_dist = glm::max(z_near_dist + 0.1f, model_furthest_dist + 1.0f); // Extend slightly behind
            }
            return glm::perspective(FOV_RAD, aspect_ratio, z_near_dist, z_far_dist);
        } else { // Orthographic projection
            float view_height = m_ortho_view_size;
            float view_width = m_ortho_view_size * aspect_ratio;

            if (m_dynamic_clipping_enabled) {
                const float FIXED_CLIP_DISTANCE = 60.0f; // 60mm from camera
                z_near_dist = FIXED_CLIP_DISTANCE;

                // Ensure z_far_dist covers the model from this new z_near_dist
                // The furthest point of the model is at -min_z_view_space distance from the camera.
                // So, the depth of the model visible *beyond* the fixed near plane is (-min_z_view_space - FIXED_CLIP_DISTANCE).
                // z_far_dist should be FIXED_CLIP_DISTANCE + (this visible depth) + padding.
                float model_depth_behind_clip_plane = glm::max(0.0f, model_furthest_dist - FIXED_CLIP_DISTANCE);
                z_far_dist = z_near_dist + model_depth_behind_clip_plane + 1.0f; // Add a small buffer
                z_far_dist = glm::max(z_near_dist + 0.1f, z_far_dist); // Ensure far is always > near
            } else {
                // Non-dynamic orthographic clipping: encompass the whole model
                z_near_dist = model_closest_dist - 1.0f; // Extend slightly in front
                z_far_dist = model_furthest_dist + 1.0f; // Extend slightly behind
                z_near_dist = glm::max(0.001f, z_near_dist); // Ensure positive near, changed from 0.01f
                z_far_dist = glm::max(z_near_dist + 0.1f, z_far_dist); // Ensure far > near
            }
            return glm::ortho(-view_width / 2.0f, view_width / 2.0f,
                               -view_height / 2.0f, view_height / 2.0f,
                               z_near_dist, z_far_dist);
        }
    }

    // Moller-Trumbore algorithm for ray-triangle intersection
    bool ray_triangle_intersect(const glm::vec3& ray_origin, const glm::vec3& ray_direction,
                                const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
                                float& out_t) {
        const float EPSILON = 0.0000001f;
        glm::vec3 edge1, edge2, h, s, q;
        float a, f, u, v;
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        h = glm::cross(ray_direction, edge2);
        a = glm::dot(edge1, h);
        if (a > -EPSILON && a < EPSILON)
            return false;    // This ray is parallel to this triangle.
        f = 1.0f / a;
        s = ray_origin - v0;
        u = f * glm::dot(s, h);
        if (u < 0.0f || u > 1.0f)
            return false;
        q = glm::cross(s, edge1);
        v = f * glm::dot(ray_direction, q);
        if (v < 0.0f || u + v > 1.0f)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        out_t = f * glm::dot(edge2, q);
        if (out_t > EPSILON) // ray intersection
            return true;
        else // This means that there is a line intersection but not a ray intersection.
            return false;
    }

    // Find the closest point on a ray to a given point
    glm::vec3 closest_point_on_ray(const glm::vec3& ray_origin, const glm::vec3& ray_direction, const glm::vec3& point) {
        glm::vec3 v = point - ray_origin;
        float t = glm::dot(v, ray_direction);
        return ray_origin + t * ray_direction;
    }

    // Sets the rotation pivot based on mouse click
    void set_rotation_pivot_from_mouse(double mouse_x, double mouse_y) {
        // Get current view and projection matrices
        glm::mat4 projection = compute_projection_matrix();
        glm::mat4 view = glm::lookAt(m_camera_position, m_camera_target, glm::vec3(0.0f, 1.0f, 0.0f));

        // Convert mouse coordinates to NDC
        // Corrected: No inversion for ndc_x here.
        float ndc_x = (2.0f * static_cast<float>(mouse_x)) / get_width() - 1.0f;
        float ndc_y = 1.0f - (2.0f * static_cast<float>(mouse_y)) / get_height(); // Y is inverted for OpenGL

        // Create ray in clip space (near plane)
        glm::vec4 ray_clip = glm::vec4(ndc_x, ndc_y, -1.0f, 1.0f);

        // Unproject to eye space
        glm::vec4 ray_eye = glm::inverse(projection) * ray_clip;
        ray_eye = glm::vec4(ray_eye.x, ray_eye.y, -1.0f, 0.0f); // Set Z to -1 for direction, W to 0 for vector

        // Unproject to world space
        glm::vec3 ray_world_dir = glm::normalize(glm::vec3(glm::inverse(view) * ray_eye));
        glm::vec3 ray_world_origin = m_camera_position;

        // Now, find intersection with geometry
        glm::vec3 intersection_point = m_model_center; // Default to center if no intersection
        float min_t = std::numeric_limits<float>::max();
        bool intersected = false;

        // Define the full model transformation matrix for consistency with rendering
        // This should match the model matrix calculation in compute_mvp()
        glm::mat4 full_model_transform = glm::translate(glm::mat4(1.0f), m_rotation_pivot) *
                                         glm::mat4_cast(m_model_rotation) *
                                         glm::translate(glm::mat4(1.0f), -m_rotation_pivot);

        // Iterate through all triangles in m_vertex_data
        // Assuming m_vertex_data stores (x,y,z,w,r,g,b,a) and m_face_index_data stores indices
        for (size_t i = 0; i < m_face_index_data.size(); i += 3) {
            GLuint idx0 = m_face_index_data[i];
            GLuint idx1 = m_face_index_data[i+1];
            GLuint idx2 = m_face_index_data[i+2];

            // Get vertices for the triangle (only position part)
            glm::vec3 v0_local(m_vertex_data[idx0 * 8], m_vertex_data[idx0 * 8 + 1], m_vertex_data[idx0 * 8 + 2]);
            glm::vec3 v1_local(m_vertex_data[idx1 * 8], m_vertex_data[idx1 * 8 + 1], m_vertex_data[idx1 * 8 + 2]);
            glm::vec3 v2_local(m_vertex_data[idx2 * 8], m_vertex_data[idx2 * 8 + 1], m_vertex_data[idx2 * 8 + 2]);

            // Transform triangle vertices by the full model transformation
            glm::vec3 v0_world = glm::vec3(full_model_transform * glm::vec4(v0_local, 1.0f));
            glm::vec3 v1_world = glm::vec3(full_model_transform * glm::vec4(v1_local, 1.0f));
            glm::vec3 v2_world = glm::vec3(full_model_transform * glm::vec4(v2_local, 1.0f));

            float t;
            if (ray_triangle_intersect(ray_world_origin, ray_world_dir, v0_world, v1_world, v2_world, t)) {
                if (t < min_t) {
                    min_t = t;
                    intersection_point = ray_world_origin + ray_world_dir * t;
                    intersected = true;
                }
            }
        }

        if (intersected) {
            m_rotation_pivot = intersection_point;
            std::cerr << "Rotation pivot set to intersection point: (" << m_rotation_pivot.x << ", " << m_rotation_pivot.y << ", " << m_rotation_pivot.z << ")" << std::endl;
        } else {
            // If no intersection, find closest vertex to ray and use that as pivot
            float min_dist_sq = std::numeric_limits<float>::max();
            glm::vec3 closest_vertex = m_model_center;
            bool found_closest_vertex = false;

            for (size_t i = 0; i < m_vertex_data.size(); i += 8) {
                glm::vec3 vertex_pos_local(m_vertex_data[i], m_vertex_data[i+1], m_vertex_data[i+2]);
                // Transform vertex by the full model transformation
                glm::vec3 vertex_pos_world = glm::vec3(full_model_transform * glm::vec4(vertex_pos_local, 1.0f));

                glm::vec3 point_on_ray = closest_point_on_ray(ray_world_origin, ray_world_dir, vertex_pos_world);
                float dist_sq = glm::distance2(vertex_pos_world, point_on_ray);

                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    closest_vertex = vertex_pos_world;
                    found_closest_vertex = true;
                }
            }

            if (found_closest_vertex) {
                m_rotation_pivot = closest_vertex;
                std::cerr << "Rotation pivot set to closest vertex: (" << m_rotation_pivot.x << ", " << m_rotation_pivot.y << ", " << m_rotation_pivot.y << ")" << std::endl;
            } else {
                // Fallback: if no geometry, pivot remains model center.
                m_rotation_pivot = m_model_center;
                std::cerr << "No geometry found for pivot, resetting to model center: (" << m_rotation_pivot.x << ", " << m_rotation_pivot.y << ", " << m_rotation_pivot.z << ")" << std::endl;
            }
        }
        queue_draw();
    }

    // Calculates the Model-View-Projection matrix
    glm::mat4 compute_mvp()
    {
        // 1. Translate model so pivot is at origin
        glm::mat4 model = glm::translate(glm::mat4(1.0f), -m_rotation_pivot);
        // 2. Apply rotation
        model = glm::mat4_cast(m_model_rotation) * model;
        // 3. Translate model back by pivot
        model = glm::translate(glm::mat4(1.0f), m_rotation_pivot) * model;

        glm::mat4 view = glm::lookAt(
            m_camera_position,
            m_camera_target,
            glm::vec3(0.0f, 1.0f, 0.0f) // World up vector
        );

        glm::mat4 projection = compute_projection_matrix(); // Use the new helper function

        return projection * view * model;
    }

    // Renders the object with hidden line removal
    void draw_object()
    {
        if (m_program == 0) return;
        if (m_vertex_data.empty() && m_edge_vertex_data.empty()) return; // Nothing to draw

        glm::mat4 mvp = compute_mvp();
        glUseProgram(m_program);
        glUniformMatrix4fv(m_mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

        // --- First pass: Render filled faces ---
        if (!m_face_index_data.empty()) {
            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE); // Enable writing to depth buffer for faces

            // Apply culling based on m_culling_enabled
            if (m_culling_enabled) {
                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);
            } else {
                glDisable(GL_CULL_FACE);
            }

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Explicitly set to fill
            glDepthFunc(GL_LESS); // Faces are rendered if they are closer

            // Disable polygon offset for faces in this pass
            glDisable(GL_POLYGON_OFFSET_FILL);
            glDisable(GL_POLYGON_OFFSET_LINE);

            glBindBuffer(GL_ARRAY_BUFFER, m_position_buffer);
            glEnableVertexAttribArray(0); // Position attribute
            glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 8, (void*)0);
            glEnableVertexAttribArray(1); // Color attribute
            glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 8, (void*)(sizeof(GLfloat) * 4));

            glUniform1i(m_is_line_pass_location, 0); // Tell shader it's not a line pass
            glUniform1i(m_apply_normal_correction_tint_location, m_apply_normal_correction_tint ? 1 : 0); // Set tint uniform
            glUniform4fv(m_normal_correction_tint_color_location, 1, glm::value_ptr(m_normal_correction_tint_color)); // Set tint color

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
            glDrawElements(GL_TRIANGLES, m_face_index_data.size(), GL_UNSIGNED_INT, 0);
        }

        // --- Second pass: Render edges ---
        if (!m_edge_line_indices.empty()) {
            glDepthMask(GL_FALSE); // DO NOT WRITE to depth buffer for lines
            glDisable(GL_CULL_FACE); // Culling does not apply to lines

            // Enable polygon offset for lines to slightly shift them towards the camera
            glEnable(GL_POLYGON_OFFSET_LINE);
            glPolygonOffset(-1.0f, -2.0f); // Adjusted offset for better line visibility

            glDepthFunc(GL_LEQUAL); // Lines are rendered if they are closer or at the same depth

            // Bind edge vertex and index buffers
            glBindBuffer(GL_ARRAY_BUFFER, m_edge_position_buffer);
            glEnableVertexAttribArray(0); // Position attribute for edges
            glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4, (void*)0);
            glDisableVertexAttribArray(1); // Edges use a uniform color, so disable vertex color attribute

            glUniform1i(m_is_line_pass_location, 1); // Tell shader it's a line pass
            glUniform4f(m_line_color_uniform_location, 0.0f, 0.0f, 0.0f, 1.0f); // Set line color to black
            glUniform1i(m_apply_normal_correction_tint_location, 0); // Disable tint for lines
            glUniform4fv(m_normal_correction_tint_color_location, 1, glm::value_ptr(m_normal_correction_tint_color)); // Still set, but not used

            glLineWidth(1.0f); // Set line width to 1 pixel

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_edge_index_buffer);
            glDrawElements(GL_LINES, m_edge_line_indices.size(), GL_UNSIGNED_INT, 0);
        }

        // --- Third pass: Render highlighted element ---
        if (m_highlighted_element.type != NONE) {
            glDisable(GL_DEPTH_TEST); // Render highlight on top
            glDisable(GL_CULL_FACE); // No culling for highlight

            glUseProgram(m_program);
            glUniformMatrix4fv(m_mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));
            glUniform1i(m_is_line_pass_location, 1); // Use line pass logic for simpler color override
            glUniform1i(m_apply_normal_correction_tint_location, 0); // Disable tint for highlights

            // Highlight color (orange)
            const float ORANGE_R = 1.0f;
            const float ORANGE_G = 0.647f;
            const float ORANGE_B = 0.0f;
            const float ORANGE_A = 1.0f;

            if (m_highlighted_element.type == FACE) {
                // Highlight faces: Render with solid color and slight offset
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(-2.0f, -4.0f); // More aggressive offset for highlight
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

                // Check the type of the highlighted face to apply different highlight color
                GeomAbs_SurfaceType highlighted_surface_type = GeomAbs_OtherSurface;
                if (m_highlighted_element.id != -1 && m_highlighted_element.id < m_face_surface_types.size()) {
                    highlighted_surface_type = m_face_surface_types[m_highlighted_element.id];
                }

                // Use a distinct color for curved surfaces (often used for fillets)
                if (highlighted_surface_type == GeomAbs_Cylinder ||
                    highlighted_surface_type == GeomAbs_Cone ||
                    highlighted_surface_type == GeomAbs_Sphere ||
                    highlighted_surface_type == GeomAbs_Torus) {
                    glUniform4f(m_line_color_uniform_location, 1.0f, 0.0f, 1.0f, 1.0f); // Magenta for curved surfaces
                } else {
                    glUniform4f(m_line_color_uniform_location, ORANGE_R, ORANGE_G, ORANGE_B, ORANGE_A); // Orange for other surfaces
                }

                glBindBuffer(GL_ARRAY_BUFFER, m_position_buffer);
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 8, (void*)0);
                glDisableVertexAttribArray(1); // Do not use vertex colors for highlight

                if (m_highlighted_element.id != -1 && m_highlighted_element.id < m_face_draw_info.size()) {
                    const auto& draw_info = m_face_draw_info[m_highlighted_element.id];
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
                    glDrawElements(GL_TRIANGLES, draw_info.second, GL_UNSIGNED_INT, (void*)(draw_info.first * sizeof(GLuint)));
                }

                glDisable(GL_POLYGON_OFFSET_FILL);

            } else if (m_highlighted_element.type == EDGE) {
                glUniform4f(m_line_color_uniform_location, ORANGE_R, ORANGE_G, ORANGE_B, ORANGE_A); // Orange highlight
                glLineWidth(3.0f); // Thicker line

                glBindBuffer(GL_ARRAY_BUFFER, m_edge_position_buffer);
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4, (void*)0);
                glDisableVertexAttribArray(1);

                if (m_highlighted_element.id != -1 && m_highlighted_element.id < m_edge_draw_info.size()) {
                    const auto& draw_info = m_edge_draw_info[m_highlighted_element.id];
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_edge_index_buffer);
                    glDrawElements(GL_LINES, draw_info.second, GL_UNSIGNED_INT, (void*)(draw_info.first * sizeof(GLuint)));
                }

            } else if (m_highlighted_element.type == VERTEX) {
                glUniform4f(m_line_color_uniform_location, ORANGE_R, ORANGE_G, ORANGE_B, ORANGE_A); // Orange highlight
                glPointSize(8.0f); // Larger point size
                glEnable(GL_PROGRAM_POINT_SIZE); // Enable point size from shader (or fixed)

                glBindBuffer(GL_ARRAY_BUFFER, m_position_buffer); // Use face vertex data for positions
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 8, (void*)0);
                glDisableVertexAttribArray(1);

                if (m_highlighted_element.id != -1 && m_highlighted_element.id * 8 < m_vertex_data.size()) {
                    glDrawArrays(GL_POINTS, m_highlighted_element.id, 1);
                }
                glDisable(GL_PROGRAM_POINT_SIZE);
            }

            glEnable(GL_DEPTH_TEST); // Re-enable depth test for next frame
        }

        // --- Reset OpenGL states to default for the next frame ---
        glDisable(GL_POLYGON_OFFSET_LINE); // Disable polygon offset for line mode
        glDepthMask(GL_TRUE); // Re-enable writing to depth buffer
        // Reset culling state to match m_culling_enabled
        if (m_culling_enabled) {
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
        } else {
            glDisable(GL_CULL_FACE);
        }
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Reset polygon mode to fill
        glLineWidth(1.0f); // Reset line width (for safety, though not critical for faces)
        glDepthFunc(GL_LESS); // Reset depth function to default GL_LESS

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glUseProgram(0);
    }
};

// Definition of generate_cube_with_opencascade
void MyGLArea::generate_cube_with_opencascade()
{
    // Clear previous data
    m_vertex_data.clear();
    m_face_index_data.clear();
    m_edge_vertex_data.clear();
    m_edge_line_indices.clear();
    m_face_draw_info.clear(); // Clear draw info
    m_edge_draw_info.clear(); // Clear draw info
    m_face_surface_types.clear(); // Clear surface types

    // Define cube vertices and colors
    // Position (x, y, z, w), Color (r, g, b, a)
    const float cube_size = 50.0f; // New variable for cube size
    const float half_size = cube_size / 2.0f;

    GLfloat vertices[] = {
        // Front face (Green)
        -half_size, -half_size,  half_size, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f, // 0
         half_size, -half_size,  half_size, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f, // 1
         half_size,  half_size,  half_size, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f, // 2
        -half_size,  half_size,  half_size, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f, // 3

        // Back face (Blue)
        -half_size, -half_size, -half_size, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f, // 4
         half_size, -half_size, -half_size, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f, // 5
         half_size,  half_size, -half_size, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f, // 6
        -half_size,  half_size, -half_size, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f, // 7

        // Top face (Red)
        -half_size,  half_size,  half_size, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f, // 8 (same as 3)
         half_size,  half_size,  half_size, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f, // 9 (same as 2)
         half_size,  half_size, -half_size, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f, // 10 (same as 6)
        -half_size,  half_size, -half_size, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f, // 11 (same as 7)

        // Bottom face (Yellow)
        -half_size, -half_size,  half_size, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f, // 12 (same as 0)
         half_size, -half_size,  half_size, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f, // 13 (same as 1)
         half_size, -half_size, -half_size, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f, // 14 (same as 5)
        -half_size, -half_size, -half_size, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f, // 15 (same as 4)

        // Right face (Cyan)
         half_size, -half_size,  half_size, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f, // 16 (same as 1)
         half_size, -half_size, -half_size, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f, // 17 (same as 5)
         half_size,  half_size, -half_size, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f, // 18 (same as 6)
         half_size,  half_size,  half_size, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f, // 19 (same as 2)

        // Left face (Magenta)
        -half_size, -half_size,  half_size, 1.0f,   1.0f, 0.0f, 1.0f, 1.0f, // 20 (same as 0)
        -half_size, -half_size, -half_size, 1.0f,   1.0f, 0.0f, 1.0f, 1.0f, // 21 (same as 4)
        -half_size,  half_size, -half_size, 1.0f,   1.0f, 0.0f, 1.0f, 1.0f, // 22 (same as 7)
        -half_size,  half_size,  half_size, 1.0f,   1.0f, 0.0f, 1.0f, 1.0f  // 23 (same as 3)
    };

    // Indices for faces (two triangles per face)
    GLuint face_indices[] = {
        // Front face
        0, 1, 2,
        0, 2, 3,

        // Back face
        4, 7, 6,
        4, 6, 5,

        // Top face
        8, 9, 10,
        8, 10, 11,

        // Bottom face
        12, 14, 13,
        12, 15, 14,

        // Right face
        16, 17, 18,
        16, 18, 19,

        // Left face
        20, 22, 21,
        20, 23, 22
    };

    // Edge vertices (just positions, will be drawn with a uniform color)
    GLfloat edge_vertices[] = {
        // Front square
        -half_size, -half_size,  half_size, 1.0f,
         half_size, -half_size,  half_size, 1.0f,
         half_size,  half_size,  half_size, 1.0f,
        -half_size,  half_size,  half_size, 1.0f,

        // Back square
        -half_size, -half_size, -half_size, 1.0f,
         half_size, -half_size, -half_size, 1.0f,
         half_size,  half_size, -half_size, 1.0f,
        -half_size,  half_size, -half_size, 1.0f,
    };

    // Indices for edges (lines)
    GLuint edge_line_indices[] = {
        // Front square
        0, 1,
        1, 2,
        2, 3,
        3, 0,

        // Back square
        4, 5,
        5, 6,
        6, 7,
        7, 4,

        // Connecting lines
        0, 4,
        1, 5,
        2, 6,
        3, 7
    };

    // Copy data to member vectors
    m_vertex_data.assign(vertices, vertices + sizeof(vertices) / sizeof(GLfloat));
    m_face_index_data.assign(face_indices, face_indices + sizeof(face_indices) / sizeof(GLuint));
    m_edge_vertex_data.assign(edge_vertices, edge_vertices + sizeof(edge_vertices) / sizeof(GLfloat));
    m_edge_line_indices.assign(edge_line_indices, edge_line_indices + sizeof(edge_line_indices) / sizeof(GLuint));

    // For cube faces:
    size_t current_face_start_index = 0;
    for (int i = 0; i < 6; ++i) { // 6 faces, 2 triangles per face, 3 indices per triangle
        m_face_draw_info.emplace_back(current_face_start_index, 6); // Each face has 6 indices (2 triangles)
        m_face_surface_types.push_back(GeomAbs_Plane); // Cube faces are planar
        current_face_start_index += 6;
    }

    // For cube edges:
    size_t current_edge_start_index = 0;
    for (int i = 0; i < 12; ++i) { // 12 edges, 2 indices per edge
        m_edge_draw_info.emplace_back(current_edge_start_index, 2);
        current_edge_start_index += 2;
    }

    // For a unit cube centered at origin, the bounding box is straightforward
    m_min_bounds = glm::vec3(-half_size, -half_size, -half_size);
    m_max_bounds = glm::vec3(half_size, half_size, half_size);
    m_model_center = (m_min_bounds + m_max_bounds) / 2.0f; // Ensure center is calculated
    m_model_dimensions = m_max_bounds - m_min_bounds;
    m_has_model_data = true;
}


class MyWindow : public Gtk::Window
{
public:
    MyWindow()
    {
        set_title("GTKmm4 OpenGL Cube/STEP Viewer - OpenCASCADE Integration with Mouse Control");
        set_default_size(800, 600); // Increase default window size

        auto main_vbox = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::VERTICAL, 12);
        set_child(*main_vbox);
        main_vbox->set_margin(12);

        m_gl_area = Gtk::make_managed<MyGLArea>();
        m_gl_area->set_hexpand(true);
        m_gl_area->set_vexpand(true);
        main_vbox->append(*m_gl_area);

        auto button_box = Gtk::make_managed<Gtk::Box>(Gtk::Orientation::HORIZONTAL, 6);
        button_box->set_halign(Gtk::Align::CENTER);
        main_vbox->append(*button_box);

        auto load_step_button = Gtk::make_managed<Gtk::Button>("Load STEP file");
        button_box->append(*load_step_button);
        load_step_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_load_step_button_clicked));

        auto toggle_culling_button = Gtk::make_managed<Gtk::Button>("Toggle Culling");
        button_box->append(*toggle_culling_button);
        toggle_culling_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_toggle_culling_button_clicked));

        auto toggle_clipping_button = Gtk::make_managed<Gtk::Button>("Toggle Dynamic Clipping");
        button_box->append(*toggle_clipping_button);
        toggle_clipping_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_toggle_clipping_button_clicked));

        auto toggle_projection_button = Gtk::make_managed<Gtk::Button>("Toggle Projection"); // New button
        button_box->append(*toggle_projection_button);
        toggle_projection_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_toggle_projection_button_clicked));

        auto front_view_button = Gtk::make_managed<Gtk::Button>("Front View"); // New button for front view
        button_box->append(*front_view_button);
        front_view_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_front_view_button_clicked));

        auto top_view_button = Gtk::make_managed<Gtk::Button>("Top View"); // New button for top view
        button_box->append(*top_view_button);
        top_view_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_top_view_button_clicked));

        auto right_view_button = Gtk::make_managed<Gtk::Button>("Right View"); // New button for right view
        button_box->append(*right_view_button);
        right_view_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_right_view_button_clicked));

        auto center_view_button = Gtk::make_managed<Gtk::Button>("Center View");
        button_box->append(*center_view_button);
        center_view_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_center_view_button_clicked));

        // New button for normal correction tint
        auto toggle_normal_tint_button = Gtk::make_managed<Gtk::Button>("Toggle Normal Tint");
        button_box->append(*toggle_normal_tint_button);
        toggle_normal_tint_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_toggle_normal_tint_button_clicked));

        auto quit_button = Gtk::make_managed<Gtk::Button>("Quit");
        button_box->append(*quit_button);
        quit_button->signal_clicked().connect(sigc::mem_fun(*this, &MyWindow::on_quit_button_clicked));
    }

protected:
    void on_load_step_button_clicked()
    {
        auto dialog = Gtk::FileDialog::create();

        // Add filters, so that only certain file types can be selected:
        auto filters = Gio::ListStore<Gtk::FileFilter>::create();

        auto filter_step = Gtk::FileFilter::create();
        filter_step->set_name("STEP files (*.step *.stp)");
        filter_step->add_pattern("*.step");
        filter_step->add_pattern("*.stp");
        filters->append(filter_step);

        auto filter_all = Gtk::FileFilter::create();
        filter_all->set_name("All  files");
        filter_all->add_pattern("*");
        filters->append(filter_all);

        dialog->set_filters(filters);

        // Show the dialog and wait for a user response asynchronously:
        // Use sigc::bind to pass the dialog itself to the finish callback
        dialog->open(sigc::bind(sigc::mem_fun(
            *this, &MyWindow::on_file_dialog_finish_step), dialog));
    }

    void on_file_dialog_finish_step(const Glib::RefPtr<Gio::AsyncResult>& result,
                                    const Glib::RefPtr<Gtk::FileDialog>& dialog)
    {
        // Handle the response:
        try
        {
            auto file = dialog->open_finish(result);

            if (file) { // Check if a file was actually selected
                // Notice that this is a std::string, not a Glib::ustring.
                auto filepath = file->get_path();
                std::cout << "STEP file selected: " <<  filepath << std::endl;
                m_gl_area->load_step_file(filepath);
            } else {
                std::cout << "No STEP file selected." << std::endl;
            }
        }
        catch (const Gtk::DialogError& err)
        {
            // Can be thrown by dialog->open_finish(result).
            std::cout << "STEP file dialog error: " << err.what() << std::endl;
        }
        catch (const Glib::Error& err)
        {
            std::cout << "Unexpected exception in STEP file dialog: " << err.what() << std::endl;
        }
    }

    void on_toggle_culling_button_clicked()
    {
        m_gl_area->toggle_culling();
    }

    void on_toggle_clipping_button_clicked() // New handler for clipping button
    {
        m_gl_area->toggle_dynamic_clipping();
    }

    void on_toggle_projection_button_clicked() // New handler for projection button
    {
        m_gl_area->toggle_projection_mode();
    }

    void on_front_view_button_clicked() // New handler for front view button
    {
        m_gl_area->set_front_view();
    }

    void on_top_view_button_clicked() // New handler for top view button
    {
        m_gl_area->set_top_view();
    }

    void on_right_view_button_clicked() // New handler for right view button
    {
        m_gl_area->set_right_view();
    }

    void on_center_view_button_clicked()
    {
        m_gl_area->reset_camera_view();
    }

    void on_toggle_normal_tint_button_clicked() // New handler for normal correction tint button
    {
        m_gl_area->toggle_normal_correction_tint();
    }

    void on_quit_button_clicked()
    {
        hide();
    }

private:
    MyGLArea* m_gl_area;
};

int main(int argc, char* argv[])
{
    auto app = Gtk::Application::create("org.gtkmm.example.OpenGLCubeStepViewer");
    return app->make_window_and_run<MyWindow>(argc, argv);
}
