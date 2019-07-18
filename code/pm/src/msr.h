
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
//* msr:                    m e a s u r e m e n t              *
//*============================================================*

#ifndef _MEASUREMENT_PM_H_
#define _MEASUREMENT_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;
class PmGeometryVector;

enum PmMeasurementType {
  PM_MEASUREMENT_UNKNOWN,
  PM_MEASUREMENT_ANGLE,
  PM_MEASUREMENT_DISTANCE,
  PM_MEASUREMENT_POSITION
  };


// PmMeasurement
// -------------

class PM_EXPORT PmMeasurement {
   public:
      PmMeasurement();
      ~PmMeasurement(){};
      void getName(string& name) { name = this->name; }
      PmMeasurementType getType() { return type; }
      static void convMeasurementType(const string str, PmMeasurementType& type);
      static PmMeasurement* create(const string name, const PmMeasurementType type);
      void setColor(const PmVector3& color);
      void getGlobalFrames(vector<bool>& frames);
      void setGlobalFrames(const vector<bool> frames);
      void getObjs(vector<void*>& objs);
      void setObjs(const vector<void*>objs);
      void getPoints(vector<PmVector3>& pts);
      void setPoints(const vector<PmVector3> pts);
      void getPrintFlag(bool& flag);
      void setPrintFlag(const bool flag);
      void setWriteResultsParams(const string file_name, const bool binary);
      virtual void writeResults(const vector<float>& time, const int last_written, 
                                const int num) = 0;
      virtual void updatePoints(const vector<PmVector3> pts)=0;
      virtual void display() = 0;
   protected:
      string name;
      PmMeasurementType type;
      vector<PmVector3> points;
      vector<PmVector3> current_points;
      vector<bool>global_frames;
      vector<void*> sobjs;
      PmGraphicsGeometry *geometry;
      PmVector3 color;

      // results i/o //
      string write_file_name;
      FILE *file_ptr;
      bool binary_format;
      bool write_init;
      int last_written;
      bool print;
   };


// PmDistMeasurement   
// ------------------------
// distance measurement     

class PM_EXPORT PmDistMeasurement : public PmMeasurement {

   public:
      PmDistMeasurement(const string name);
      ~PmDistMeasurement(){};
      void display();
      void updatePoints(const vector<PmVector3> pts);
      void writeResults(const vector<float>& time, const int last_written, 
                        const int num);
   private:
     vector<float> distances;
   };


// PmAngleMeasurement
// ------------------------
// distance measurement

class PM_EXPORT PmAngleMeasurement : public PmMeasurement {

   public:
      PmAngleMeasurement(const string name);
      ~PmAngleMeasurement(){};
      void display();
      void updatePoints(const vector<PmVector3> pts);
      void writeResults(const vector<float>& time, const int last_written, 
                        const int num);
      void setVectors(PmGeometryVector *vgeom1, PmBody *body1, PmGeometryVector *vgeom2, 
                      PmBody *body2);
   private:
     vector<float> angles;
     PmGeometryVector *vector_geom1, *vector_geom2;
     PmBody *vector_body1, *vector_body2;
   };

// PmPositionMeasurement
// ---------------------
// position measurement

class PM_EXPORT PmPositionMeasurement : public PmMeasurement {

   public:
      PmPositionMeasurement(const string name);
      ~PmPositionMeasurement(){};
      void display();
      void updatePoints(const vector<PmVector3> pts);
      void writeResults(const vector<float>& time, const int last_written,
                        const int num);
   private:
      vector<PmVector3> positions;
   };

}

#endif

