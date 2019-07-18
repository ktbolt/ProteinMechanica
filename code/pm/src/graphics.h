/* -------------------------------------------------------------------------- *
 *                              Protein Mechanica                             *
 * -------------------------------------------------------------------------- *
 * This is part of the Protein Mechanica coarse-grained molecular motor       *
 * modeling application originating from Simbios, the NIH National Center for *
 * Physics-Based Simulation of Biological Structures at Stanford, funded      *
 * under the NIH Roadmap for Medical Research, grant U54 GM072970.            *
 * See https://simtk.org.                                                     *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: David Parker                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//*============================================================*
//* gr:                 g r a p h i c s                        *
//*============================================================*

#ifndef _GRAPHICS_PM_H_ 
#define _GRAPHICS_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

class PmMolecule;

extern const char* pmGraphicsGeomAtoms;
extern const char* pmGraphicsGeomBackbone;
extern const char* pmGraphicsGeomBackboneTube;
extern const char* pmGraphicsGeomBackbonePlanes;
extern const char* pmGraphicsGeomBonds;
extern const char* pmGraphicsGeomHbonds;

typedef enum {
  PM_GEOMETRY_DISPLAY_UNKNOWN,
  PM_GEOMETRY_DISPLAY_POINT,
  PM_GEOMETRY_DISPLAY_LINE,
  PM_GEOMETRY_DISPLAY_SOLID,
  PM_GEOMETRY_DISPLAY_OUTLINE,
  PM_GEOMETRY_DISPLAY_SIZE 
  } PmGeometryDisplayType;

typedef enum {
  PM_GEOMETRY_SHADING_UNKNOWN,
  PM_GEOMETRY_SHADING_NONE,
  PM_GEOMETRY_SHADING_FLAT,
  PM_GEOMETRY_SHADING_COLOR,
  PM_GEOMETRY_SHADING_NORMAL,
  PM_GEOMETRY_SHADING_SIZE  
  } PmGeometryShadingType;

typedef struct PmEllipsoidData {
  float radius1, radius2, radius3;
  PmVector3 axis1, axis2, axis3;
  } PmEllipsoidData;

// PmGraphicsPickPoints
// ----------------------

class PmGraphicsPickPoints {
  public:
     PmGraphicsPickPoints();
     void requestPoints(const int num_pts, const string cmd);
     void addPoint(const PmVector3& pt);
     bool isActive();

  private:
    bool active;
    string command;
    int num_points;
    int num_request_pts;
    vector<PmVector3> points;
  };


// PmGraphicsAttributes
// --------------------

class PmGraphicsAttributes {
  public:
    bool color_set; PmVector3 color;
    bool line_width_set; float line_width;
    bool disjoint_set; bool disjoint; 
    bool display_type_set; PmGeometryDisplayType display_type;
    bool marker_set; bool marker;
    bool scale_set; float scale;
    bool shading_type_set; PmGeometryShadingType shading_type;
    bool visible_set; bool visible; 
    bool visible_aux_set; bool visible_aux; 
    bool lighting_set, twosided_lighting;

    PmGraphicsAttributes() { init(); }

    PmGraphicsAttributes& operator=(const PmGraphicsAttributes& src) {
      if (src.color_set) color = src.color;
      if (src.disjoint_set) disjoint = src.disjoint; 
      if (src.display_type_set) display_type = src.display_type;
      if (src.line_width_set) line_width = src.line_width;
      if (src.marker_set) marker = src.marker;
      if (src.scale_set) scale = src.scale;
      if (src.shading_type_set) shading_type = src.shading_type;
      if (src.visible_set) visible = src.visible;
      if (src.visible_aux_set) visible_aux = src.visible_aux;
      if (src.lighting_set) twosided_lighting = src.twosided_lighting;
      return *this;
      }

    void setColor(const PmVector3& color) {
      color_set = true;
      this->color = color;
      }

    void setDisjoint(const bool flag) {
      disjoint_set = true;
      this->disjoint = flag;
      }

    void setDisplayType(const PmGeometryDisplayType type) {
      display_type_set = true;
      this->display_type = type;
      }

    void setLighting(const bool flag) {
      lighting_set = true;
      this->twosided_lighting = flag;
      }

    void setLineWidth(const float width) {
      line_width_set = true;
      this->line_width = width;
      }

    void setMarker(const bool marker) {
      marker_set = true;
      this->marker = marker;
      }

    void setScale(const float scale) {
      scale_set = true;
      this->scale = scale;
      }

    void setShadingType(const PmGeometryShadingType type) {
      shading_type_set = true;
      this->shading_type = type;
      }

    void setVisible(const bool visible) {
      visible_set = true;
      this->visible = visible;
      }

    void setVisibleAux(const bool val) {
      visible_aux_set = true;
      this->visible_aux = val;
      }

    void init() {
      color_set = false;
      color.set(1,1,1);
      line_width_set = false;
      line_width = 1.0;
      scale_set = false;
      scale = 1.0;
      disjoint = false;
      display_type_set = false;
      display_type = PM_GEOMETRY_DISPLAY_SOLID;
      shading_type_set = false;
      shading_type = PM_GEOMETRY_SHADING_COLOR;
      visible_set = true;
      visible = true;
      visible_aux_set = false;
      visible_aux = false;
      lighting_set = false;
      twosided_lighting = false;
      }
  };


// PmGrInterface 
// -------------

 class PM_EXPORT PmGraphicsInterface {

    public:

      PmGraphicsInterface (bool graphics_enabled);
      ~PmGraphicsInterface ();
      void getContext (void **ptr);
      void init();
      void setCenter(PmVector3& center);
      void procEvents(const char *script);
      void setExtent(PmExtent& extent);
      void initPicking();
      void displayPick(PmVector3& point);
      void getCurrentPick(PmVector3& point);
      void rotateScene(PmVector3& rot);
      void scaleScene(const float scale);
      void translateScene(PmVector3& vec);
      void getSceneXform(PmXform& xform);
      void update();
      void setWindowColor(const PmVector3& color);
      void recordWindow(const string fname, const string format);
      void writeWindow(const string fname, const string format);
      void requestPickPoints(const int num_pts, const string cmd);

      static PmGeometryDisplayType mapRenderType(PmMoleculeRenderType rtype);

      static void mapDataToColors(float vmin, float vmax, int num_data, float *data, 
                                  PmVector3 *colors);
      static void getSpectrumColorMap(int num_colors, vector<PmVector3>& colors);

    private:

      int initialized;
      bool graphics_enabled;
      int win_width, win_height;
      void *prv_data;

      class PmGraphicsPoint *pick_geom;
      PmVector3 current_picked_point;
      PmGraphicsPickPoints pick_points;

      void api_init(); 
      void api_create_window (char *name, int x, int y, int win_width, int win_height);
      void api_proc_events (const char *promp, const char *script);
      void api_set_extent (PmExtent& extent);
      void api_get_context (void **ptr);
   };


// PmGraphicsGeometry 
// ------------------

 class PM_EXPORT PmGraphicsGeometry {
   public:
     string name;
     PmMolecule *mol;
     void *context;
     PmGraphicsAttributes attributes;

     bool copy_vertices;
     int num_vertices;
     PmVector3 *vertices;
     int num_colors;
     PmVector3 *colors;
     PmXform xform;

     PmGraphicsGeometry();
     PmGraphicsGeometry(const string name);

     virtual void display() = 0;
     virtual void create() = 0;
     void setAttributes(PmGraphicsAttributes& atts);
     void setColors(int num_colors, PmVector3 *colors);
     void setColors(vector<PmVector3>& colors);
     bool hasColors () { return this->colors != NULL; }
     void setXform(PmXform& xform) { this->xform = xform; }
     void update(int num_verts, PmVector3 *verts);
     void updateColors(int num_colors, PmVector3 *colors);

     static void convDisplayType (const string s, PmGeometryDisplayType& dtype);
     static void convShadingType (const string s, PmGeometryShadingType& stype);
   };
   

// PmGraphicsAtoms
// ---------------
// atom geometry

 class PM_EXPORT PmGraphicsAtoms: public PmGraphicsGeometry {
   public:
      PmGraphicsAtoms(const string name, int num_verts, PmVector3 *verts, 
                      PmVector3 *colors);
      void create(){};
      void display();
      void setDisplay(PmMoleculeDisplayType type) { disp_type = type; }
      void setRadius(int num, float *rads) { radii = rads; }
      void setRender (PmMoleculeRenderType type) { ren_type = type; }

    private:
      float *radii;
      PmMoleculeDisplayType disp_type;
      PmMoleculeRenderType ren_type;
   };


// PmGraphicsAxes
// --------------
// axes geometry

 class PM_EXPORT PmGraphicsAxes : public PmGraphicsGeometry {
   public:
     PmGraphicsAxes(const string name, PmVector3 center, PmVector3 axis1, 
                     PmVector3 axis2, PmVector3 axis3);
     void display();
     void create(){};
   };


// PmGraphicsBackbone 
// ------------------
// backbone geometry 

 class PM_EXPORT PmGraphicsBackbone: public PmGraphicsGeometry {
   public:
      PmGraphicsBackbone();
      PmGraphicsBackbone(const string name, int num_verts, PmVector3 *verts);
      PmGraphicsBackbone(const string name, int num, PmConn *conn, int num_verts, 
                         PmVector3 *verts);
      void create();
      void display();
    private:
      int number;
      PmConn *connectivity;
   };

// PmGraphicsBackboneTube
// ----------------------
// backbone tube geometry

 class PM_EXPORT PmGraphicsBackboneTube: public PmGraphicsGeometry {
   public:
      PmGraphicsBackboneTube();
      PmGraphicsBackboneTube(const string name, int num_verts, PmVector3 *verts, 
                             int res, float width);
      PmGraphicsBackboneTube(const string name, int num, vector<int>ends, int num_verts,
                             PmVector3 *verts, int res, float width);
      void create();
      void createFit(int count, int num_backbone_verts, PmVector3 *backbone_verts);
      void display();
    private:
      int number;
      int resolution;
      float width;
      PmConn *connectivity;
      vector<int> ends;
   };


// PmGraphicsBonds
// ---------------------
// atomic bonds geometry

 class PM_EXPORT PmGraphicsBonds : public PmGraphicsGeometry {
   public:
      PmGraphicsBonds(const string name, int num_verts, PmVector3 *verts, 
                      int num_bonds, PmConn *bonds, PmVector3 *colors);

      void setAtomColorType (PmAtomColorType type) { this->atom_color_type = type; }
      void setAtomColor (PmVector3& color) { this->atom_color = color; }
      void setBondAtoms (bool flag) { this->bond_atoms = flag; }
      void setBondColor (PmVector3& color) { this->bond_color = color; }
      void setBondColorType (PmAtomColorType type) { this->bond_color_type = type;}
      void setDisplayBondAtoms (bool flag) { this->bond_atoms = flag; }
      void setAtomRenderType (PmMoleculeRenderType type) { this->atom_render = type; }

      void create(){};
      void display();

    private:
      PmGraphicsAtoms *atom_geom; 
      void *sphere_context;
      int num_bonds;
      PmConn *bonds;
      PmAtomColorType atom_color_type, bond_color_type;
      PmMoleculeRenderType atom_render;
      PmVector3 atom_color, bond_color; 
      bool bond_atoms, base_planes, cbonly;
      float atom_radius;
   };

// PmGraphicsPlanes
// ---------------------
// atomic planes geometry

 class PM_EXPORT PmGraphicsPlanes: public PmGraphicsGeometry {
   public:
      PmGraphicsPlanes(const string name, int num_verts, PmVector3 *verts,
                       int num_planes, PmConn *conn, PmVector3 *colors);

      void setAtomColorType (PmAtomColorType type) { this->atom_color_type = type; }
      void setAtomColor (PmVector3& color) { this->atom_color = color; }
      void setPlaneAtoms (bool flag) { this->bond_atoms = flag; }
      void setPlaneColor (PmVector3& color) { this->bond_color = color; }
      void setPlaneColorType (PmAtomColorType type) { this->bond_color_type = type;}
      void setDisplayPlaneAtoms (bool flag) { this->bond_atoms = flag; }
      void setAtomRenderType (PmMoleculeRenderType type) { this->atom_render = type; }

      void create(){};
      void display();

    private:
      PmGraphicsAtoms *atom_geom;
      void *sphere_context;
      int num_planes;
      PmConn *connectivity;
      PmAtomColorType atom_color_type, bond_color_type;
      PmMoleculeRenderType atom_render;
      PmVector3 atom_color, bond_color;
      bool bond_atoms, base_planes, cbonly;
      float atom_radius;
   };

// PmGraphicsEllipsoid
// -------------------
// ellipsoid geometry

 class PM_EXPORT PmGraphicsEllipsoid : public PmGraphicsGeometry {
   public:
      PmGraphicsEllipsoid();
      PmGraphicsEllipsoid(const string name, int num_verts, PmVector3 *verts,
                          PmEllipsoidData *data);
                          
      void create(){};
      void display();
    private:
      PmEllipsoidData *axes_data;
   };

// PmGraphicsCylinder
// ------------------
// cylinder geometry

 class PM_EXPORT PmGraphicsCylinder : public PmGraphicsGeometry {
   public:
     PmGraphicsCylinder(const string name, PmVector3 center, float radius, float length,
                        PmVector3 axis);
     void create(){};
     void display();
   private:
     PmVector3 center;
     float radius, length; 
     PmVector3 axis;
   };


// PmGraphicsLine
// --------------
// line geometry

 class PM_EXPORT PmGraphicsLine : public PmGraphicsGeometry {
   public:
      PmGraphicsLine();
      PmGraphicsLine(const string name, int num_verts, PmVector3 *verts, bool copy=true);
      void create(){};
      void display();
      void setDisjoint(const bool flag);
   };

// PmGraphicsPoint
// ---------------
// point geometry

 class PM_EXPORT PmGraphicsPoint: public PmGraphicsGeometry {
   public:
      PmGraphicsPoint();
      PmGraphicsPoint(const string name, int num_verts, PmVector3 *verts);
      void create(){};
      void display();
    private:
   };

// PmGraphicsPolygon
// -----------------
// polygon geometry

 class PM_EXPORT PmGraphicsPolygon : public PmGraphicsGeometry {
   public:
      PmGraphicsPolygon();
      PmGraphicsPolygon(const string name, int num_poly, PmConn *conn, int hgeom, 
                        int num_verts, PmVector3 *verts);
      void create(){};
      void display();
   private:
      int num_polygons, nodes_per_poly; 
      PmConn *connectivity;
   };

// PmGraphicsRestraint
// -------------------
// restraint geometry

 class PM_EXPORT PmGraphicsRestraint : public PmGraphicsGeometry {
   public:
     PmGraphicsRestraint(const string name, PmVector3& pt1, PmVector3& pt2);
     void display();
     void create(){};
     void update(PmVector3& pt1, PmVector3& pt2);
   };

// PmGraphicsSphere
// ----------------
// sphere geometry

 class PM_EXPORT PmGraphicsSphere: public PmGraphicsGeometry {
   public:
      PmGraphicsSphere();
      PmGraphicsSphere(const string name, int num_verts, PmVector3 *verts, 
                       float *radius);
      void create(){};
      void display();
    private:
      float *radius;
   };

// PmGraphicsSurface 
// -----------------
// surface geometry

 class PM_EXPORT PmGraphicsSurface : public PmGraphicsGeometry {
   public:
      PmGraphicsSurface(const string name, int num_verts, PmVector3 *verts, 
                        int num_polys, PmConn *conn, int hgeom);
      virtual void create();
      virtual void display();
    private:
      int num_polygons;
      PmConn *connectivity;
      int hgeom;
   };


// PmGraphicsJoint 
// -------------------
// joint geometry

 class PM_EXPORT PmGraphicsJoint {
   public:
      PmGraphicsJoint() {};
      PmGraphicsJoint(const string name, int n, PmVector3 *verts);
      void setAttributes(PmGraphicsAttributes& atts) { this->attributes = atts; }
      virtual void create() = 0;
      virtual void display() = 0;
      void update(const PmVector3& pos);
    protected:
      string name;
      void *context, *aux_context;
      PmXform xform;
      int num_vertices;
      PmVector3 *vertices;
      PmGraphicsAttributes attributes;
   };


// PmGraphicsBallJoint
// -------------------
// ball joint geometry

 class PM_EXPORT PmGraphicsBallJoint : public PmGraphicsJoint {
   public:
      PmGraphicsBallJoint(const string name, PmVector3 *verts, PmVector3 a[3]) :
        PmGraphicsJoint(name, 1, verts) { this->axes[0] = a[0]; this->axes[1] = a[1]; 
                                          this->axes[2] = a[2]; }
      void create();
      void display();
      void xform(PmXform& xform);
   private:
     PmVector3 axes[3];
   };

// PmGraphicsFreeJoint
// -------------------
// free joint geometry

 class PM_EXPORT PmGraphicsFreeJoint : public PmGraphicsJoint {
   public:
      PmGraphicsFreeJoint(const string name, PmVector3 *verts) :
        PmGraphicsJoint(name, 1, verts) { };
      void create();
      void display();
   };


// PmGraphicsWeldJoint
// ---------------------
// ground joint geometry

 class PM_EXPORT PmGraphicsWeldJoint : public PmGraphicsJoint {
   public:
      PmGraphicsWeldJoint(const string name, PmVector3 *verts) :
        PmGraphicsJoint(name, 1, verts) { }
      void create();
      void display();
      void xform(PmXform& xform);
   private:
   };


// PmGraphicsHingeJoint
// -------------------
// pin joint geometry

 class PM_EXPORT PmGraphicsHingeJoint : public PmGraphicsJoint {
   public:
      PmGraphicsHingeJoint(const string name, int n, PmVector3 *verts, PmVector3& axis) :
        PmGraphicsJoint(name, n, verts) { this->axis = axis; }
      void create();
      void display();
      void xform(PmXform& xform);
   private:
     PmVector3 axis;
   };


// PmGraphicsUniversalJoint
// ------------------------
// universal joint geometry

 class PM_EXPORT PmGraphicsUniversalJoint : public PmGraphicsJoint {
   public:
      PmGraphicsUniversalJoint(const string name, PmVector3 *verts, PmVector3 a[2]) :
        PmGraphicsJoint(name, 1, verts) { this->axes[0] = a[0]; this->axes[1] = a[1]; }
      void create();
      void display();
      void xform(PmXform& xform);
   private:
     PmVector3 axes[2];
   };


// PmGraphicsModelJoint
// --------------------
// model joint geometry

 class PM_EXPORT PmGraphicsModelJoint {
   public:
      PmGraphicsModelJoint (const string name, int num_coords, PmVector3 *coords);
      void setAttributes(PmGraphicsAttributes& atts) { this->attributes = atts; }
      void create();
      void display();
   private:
      string name;
      void *context;
      int num_vertices;
      PmVector3 *vertices;
      PmGraphicsAttributes attributes;
   };


// PmGraphicsModelBody
// ---------------------
// model body geometry

 class PM_EXPORT PmGraphicsModelBody {
   public:
      PmGraphicsModelBody (const string name, int num_coords, PmVector3 *coords);
      void setAttributes(PmGraphicsAttributes& atts) { this->attributes = atts; }
      void create();
      void display();
   private:
     string name;
     void *context, *sphere_context;
     int num_vertices;
     PmVector3 *vertices;
     PmGraphicsAttributes attributes;
   };

// PmGraphicsGrid
// --------------
// grid 

 class PM_EXPORT PmGraphicsGrid {
   public:
      PmGraphicsGrid (const string name, PmVector3 center, PmVector3 number, 
                      PmVector3 spacing);
      void setAttributes(PmGraphicsAttributes& atts) { this->attributes = atts; }
      void display(const bool show);

   private:
     string name;
     void *context; 
     PmVector3 center; 
     PmVector3 number;
     PmVector3 spacing;
     int num_vertices;
     PmVector3 *vertices;
     PmGraphicsAttributes attributes;
   };

// PmGraphicsPlane
// ---------------

 class PM_EXPORT PmGraphicsPlane {
   public:
      PmGraphicsPlane (const string name, PmVector3 center, PmVector3 normal,
                       int number, float spacing);
      void setAttributes(PmGraphicsAttributes& atts) { this->attributes = atts; }
      void display(const bool show);

   private:
     string name;
     void *context;
     PmVector3 center, normal;
     int number;
     float spacing;
     int num_vertices;
     PmVector3 *vertices;
     int num_polygons;
     PmConn *connectivity;
     PmGraphicsAttributes attributes;
   };

}

#endif

