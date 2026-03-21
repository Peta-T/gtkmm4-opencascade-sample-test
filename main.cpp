// many thanks to Kirill Gavrilov for OcctGlTools and OcctGtkTools

 //   g++ -std=c++17 -g main.cpp OcctGtkGLAreaViewer.cpp OcctGtkWindowSample.cpp OcctGlTools.cpp OcctGtkTools.cpp OcctInputBridge.cpp AdvancedMouseTracker.cpp app_resource.o -o gstepview.exe -mwindows $(pkg-config --cflags --libs gtkmm-4.0 epoxy yaml-cpp librsvg-2.0) -I./ -I/ucrt64/include/opencascade -I/ucrt64/include/opencascade/Standard -D_USE_MATH_DEFINES -lopengl32 -lgdi32 -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep -lTKCAF -lTKCDF -lTKDCAF -lTKDECascade -lTKDEGLTF -lTKDEIGES -lTKDEOBJ -lTKDEPLY -lTKDE -lTKDESTEP -lTKDESTL -lTKDEVRML -lTKDraw -lTKernel -lTKExpress -lTKFeat -lTKFillet -lTKG2d -lTKG3d -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKLCAF -lTKMath -lTKMesh -lTKMeshVS -lTKOffset -lTKOpenGl -lTKOpenGlTest -lTKPrim -lTKQADraw -lTKRWMesh -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKTObjDRAW -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest -lTKXCAF -lTKXDEDRAW -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj -lTKXmlXCAF -lTKXSBase -lTKXSDRAWDE -lTKXSDRAWGLTF -lTKXSDRAWIGES -lTKXSDRAWOBJ -lTKXSDRAWPLY -lTKXSDRAW -lTKXSDRAWSTEP -lTKXSDRAWSTL -lTKXSDRAWVRML

#include "OcctGtkWindowSample.h"

#include "OcctGtkTools.h"

#include <Message.hxx>
#include <OSD.hxx>
#include <OSD_Environment.hxx>

#include <gtkmm.h>

int main(int theNbArgs, char* theArgVec[])
{
    Glib::setenv("GTK_THEME", "Adwaita", true);
    Glib::setenv("GSK_RENDERER", "cairo", true);
    putenv((char*)"LC_ALL=C");
    putenv((char*)"LANG=C");

    OSD::SetSignal(false);

    OcctGtkTools::gtkGlPlatformSetup();

    // 1. CHANGE: Add the HANDLES_OPEN flag
    Glib::RefPtr<Gtk::Application> aGtkApp =
      Gtk::Application::create("com.petat.gstepview",
                               Gio::Application::Flags::HANDLES_OPEN | Gio::Application::Flags::NON_UNIQUE);

    OcctGtkWindowSample* mainWindow = nullptr;

    // 2. Handler for standard launch (WITHOUT arguments)
    aGtkApp->signal_activate().connect([&aGtkApp, &mainWindow]() {
        if (!mainWindow) {
            mainWindow = new OcctGtkWindowSample();
            aGtkApp->add_window(*mainWindow);
            // Memory management: delete the window after it is closed
            mainWindow->signal_hide().connect([mainWindow]() { delete mainWindow; });
        }
        mainWindow->present();
    });

    // 3. Handler for launching with a file (WITH arguments)
    aGtkApp->signal_open().connect([&aGtkApp, &mainWindow](const std::vector<Glib::RefPtr<Gio::File>>& files, const Glib::ustring& /*hint*/) {
        if (!mainWindow) {
            mainWindow = new OcctGtkWindowSample();
            aGtkApp->add_window(*mainWindow);
            mainWindow->signal_hide().connect([mainWindow]() { delete mainWindow; });
        }
        mainWindow->present();

        if (!files.empty()) {
            Glib::ustring filePath = files[0]->get_path();

            // We start loading only when the application is idle.
            // This ensures that the gray UI of the application is drawn first and doesn't appear frozen at startup.
            Glib::signal_idle().connect_once([mainWindow, filePath]() {
                mainWindow->openFile(filePath);
            });
        }
    });

    return aGtkApp->run(theNbArgs, theArgVec);
}