# gtkmm4 easy Application with Opencascade sample test step file viewer
This simplified sample C++ code which try show possible using gtkmm4  with integraton Opencascade this code is possible to build and run on Debian13 or on Windows in MSYS2 enviroment.

Building on Windows:
====================

1.) Install MSYS2
-----------------
Download and run the msys2 installer from http://msys2.github.io/ Tested only with the 64bit version. Make sure that the path you select for installation doesn’t contain any spaces.

2.) Start MSYS console
----------------------
Launch the Start Menu item “MSYS2 mingw 64 bit” you should be greeted with a console window. All steps below refer to what you should type into that window.

3.) Install updates
-------------------
Type:

   pacman -Syu

if it tells you to close restart msys, close the console window and start it again. Then run pacman -Syu again.

4.) Install dependencies
------------------------
Type/paste

   pacman -S \\ \
   mingw-w64-x86_64-gcc \\ \
   mingw-w64-x86_64-pkgconf \\ \
   mingw-w64-x86_64-gtkmm4 \\ \
   zip \\ \
   unzip \\ \
   git \\ \
   mingw-w64-x86_64-opencascade \\ \
   --needed

When prompted, just hit return. Sit back and wait for it to install what’s almost a complete linux environment.

Before continuing you may change to another directory. It easiest to type cd followed by a space and drop the folder you want to change to on the window.

5.) Clone gtkmm4-opencascade-sample-test by type/paste on commandline:
---------------------------------------------------------------------

   git clone https://github.com/Peta-T/gtkmm4-opencascade-sample-test \
   cd gtkmm4-opencascade-sample-test

6.) Build it - copy/paste on command line:
------------------------------------

   g++ -std=c++17 -g main.cc -o app $(pkg-config --cflags --libs gtkmm-4.0 epoxy librsvg-2.0) -I/usr/include/gtk-4.0 -I/usr/include/opencascade -I/usr/include/opencascade/Standard -I/mingw64/include/gtk-4.0 -I/mingw64/include/opencascade -I/mingw64/include/opencascade/Standard -I/usr/include/gtk-4.0/gdk/x11 -L/usr/local/lib/x86_64-linux-gnu -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep -lTKCAF -lTKCDF -lTKDCAF -lTKDECascade -lTKDEGLTF -lTKDEIGES -lTKDEOBJ -lTKDEPLY -lTKDE -lTKDESTEP -lTKDESTL -lTKDEVRML -lTKDraw -lTKernel -lTKExpress -lTKFeat -lTKFillet -lTKG2d -lTKG3d -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKLCAF -lTKMath -lTKMesh -lTKMeshVS -lTKOffset -lTKOpenGl -lTKOpenGlTest -lTKPrim -lTKQADraw -lTKRWMesh -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKTObjDRAW -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest -lTKXCAF -lTKXDEDRAW -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj -lTKXmlXCAF -lTKXSBase -lTKXSDRAWDE -lTKXSDRAWGLTF -lTKXSDRAWIGES -lTKXSDRAWOBJ -lTKXSDRAWPLY -lTKXSDRAW -lTKXSDRAWSTEP -lTKXSDRAWSTL -lTKXSDRAWVRML  -D_USE_MATH_DEFINES

7.) Run app - type on command line:
-----------------------------------

   ./app


