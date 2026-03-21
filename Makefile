CXX = g++
CXXFLAGS = -std=c++17 -g $(shell pkg-config --cflags gtkmm-4.0 epoxy yaml-cpp librsvg-2.0) -I/usr/include/opencascade
LIBS = $(shell pkg-config --libs gtkmm-4.0 epoxy yaml-cpp librsvg-2.0) \
       -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep \
       -lTKCAF -lTKCDF -lTKDCAF -lTKDECascade -lTKDEGLTF -lTKDEIGES \
       -lTKDEOBJ -lTKDEPLY -lTKDE -lTKDESTEP -lTKDESTL -lTKDEVRML \
       -lTKDraw -lTKernel -lTKExpress -lTKFeat -lTKFillet -lTKG2d -lTKG3d \
       -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKLCAF -lTKMath -lTKMesh \
       -lTKMeshVS -lTKOffset -lTKOpenGl -lTKOpenGlTest -lTKPrim -lTKQADraw \
       -lTKRWMesh -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKTObjDRAW \
       -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest \
       -lTKXCAF -lTKXDEDRAW -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj \
       -lTKXmlXCAF -lTKXSBase -lTKXSDRAWDE -lTKXSDRAWGLTF -lTKXSDRAWIGES \
       -lTKXSDRAWOBJ -lTKXSDRAWPLY -lTKXSDRAW -lTKXSDRAWSTEP -lTKXSDRAWSTL -lTKXSDRAWVRML \
       -lGL -lEGL

SRCS = main.cpp OcctGtkGLAreaViewer.cpp OcctGtkWindowSample.cpp OcctGlTools.cpp \
       OcctGtkTools.cpp OcctInputBridge.cpp AdvancedMouseTracker.cpp

all: gstepview

gstepview: $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o gstepview $(LIBS)

install:
	# Instalace binárního souboru
	install -D gstepview $(DESTDIR)/usr/bin/gstepview
	# Instalace spouštěče do menu
	install -m 644 -D gstepview.desktop $(DESTDIR)/usr/share/applications/gstepview.desktop
	# Instalace ikony pro systém
	install -m 644 -D app_icon.png $(DESTDIR)/usr/share/pixmaps/gstepview.png

clean:
	rm -f gstepview
