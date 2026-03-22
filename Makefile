# Název výsledné aplikace
APP_NAME = gstepview

# Detekce OS (Windows_NT je standardní proměnná prostředí na Windows)
ifeq ($(OS),Windows_NT)
    # --- Nastavení pro MSYS2 / Windows ---
    EXE_EXT = .exe
    # Cesta k OpenCASCADE v MSYS2 UCRT64
    OCCT_INC = -I/ucrt64/include/opencascade
    # Knihovny specifické pro Windows
    PLATFORM_LIBS = -lopengl32 -lEGL -lgdi32
    # Flag pro skrytí konzole (volitelné, odstraňte pro ladění)
    WIN_FLAGS = -mwindows
    # Resource soubor s ikonou
    RESOURCES = app_resource.o
else
    # --- Nastavení pro LINUX ---
    EXE_EXT =
    # Na Linuxu bývá OCCT v /usr/include/opencascade
    OCCT_INC = -I/usr/include/opencascade
    PLATFORM_LIBS = -lGL -lEGL
    WIN_FLAGS =
    RESOURCES =
endif

# Kompilátor a základní příznaky
CXX = g++
CXXFLAGS = -std=c++17 -g -D_USE_MATH_DEFINES

# Získání příznaků pro GTK4 a další knihovny přes pkg-config
GTK_CFLAGS = $(shell pkg-config --cflags gtkmm-4.0 epoxy yaml-cpp librsvg-2.0)
GTK_LIBS   = $(shell pkg-config --libs gtkmm-4.0 epoxy yaml-cpp librsvg-2.0)

# Seznam všech zdrojových souborů [cite: 1]
SOURCES = main.cpp \
          OcctGtkGLAreaViewer.cpp \
          OcctGtkWindowSample.cpp \
          OcctGlTools.cpp \
          OcctGtkTools.cpp \
          OcctInputBridge.cpp \
          AdvancedMouseTracker.cpp

# Seznam všech OpenCASCADE knihoven (TK*)
OCCT_LIBS = -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep \
            -lTKCAF -lTKCDF -lTKDCAF -lTKDECascade -lTKDEGLTF -lTKDEIGES \
            -lTKDEOBJ -lTKDEPLY -lTKDE -lTKDESTEP -lTKDESTL -lTKDEVRML \
            -lTKDraw -lTKernel -lTKExpress -lTKFeat -lTKFillet -lTKG2d -lTKG3d \
            -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKLCAF -lTKMath -lTKMesh \
            -lTKMeshVS -lTKOffset -lTKOpenGl -lTKOpenGlTest -lTKPrim -lTKQADraw \
            -lTKRWMesh -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKTObjDRAW \
            -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest \
            -lTKXCAF -lTKXDEDRAW -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj \
            -lTKXmlXCAF -lTKXSBase -lTKXSDRAWDE -lTKXSDRAWGLTF -lTKXSDRAWIGES \
            -lTKXSDRAWOBJ -lTKXSDRAWPLY -lTKXSDRAW -lTKXSDRAWSTEP \
            -lTKXSDRAWSTL -lTKXSDRAWVRML

# Finální cílový soubor
TARGET = $(APP_NAME)$(EXE_EXT)

# Pravidlo pro vše (výchozí)
all: $(TARGET)

# Sestavení aplikace
$(TARGET): $(SOURCES) $(RESOURCES)
	$(CXX) $(CXXFLAGS) $(OCCT_INC) $(GTK_CFLAGS) $^ -o $@ \
	$(WIN_FLAGS) $(GTK_LIBS) $(OCCT_LIBS) $(PLATFORM_LIBS)

# Pravidlo pro Windows Resource soubor (ikona)
app_resource.o: app_resource.rc
	windres app_resource.rc -o app_resource.o

# Čištění projektu
clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
