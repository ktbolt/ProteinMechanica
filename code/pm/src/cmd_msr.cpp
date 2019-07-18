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

////////////////////////////////////////////////////////////////
//        m e a s u r e m e n t    c o m m a n d s           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========          pm_CmdMeasure                 ==========*
//*============================================================*
// process measure command.

void
pm_CmdMeasure (PmCmdDataList& dlist)
  {
  PmCmdData data, pdata1, pdata2, vdata;
  string dv, pstr1, pstr2, vstr; 
  vector<PmVector3> points;
  PmVector3 point, v, dir, vec;
  float dist;
  static PmGraphicsLine *dist_geom = NULL; 
  PmGraphicsAttributes atts;
  PmVector3 color;
  static FILE *fp=NULL;
  string file_name;
  bool pick, in_dir=false;

  dlist.getNext (data);

  //===== measure distance =====//

  if (data.name == "distance") {
    fprintf (stderr, "\n>>> distance measurement \n"); 

    if (dlist.getString("file", file_name)) {
      bool bval;

      if (dlist.getBoolean("file_initialize", bval)) {
        fp = fopen (file_name.c_str(), "w");
        fclose(fp);
        fp = fopen (file_name.c_str(), "a");
        return;
        }
      }
    else {
      fp = NULL;
      }

    if (dlist.getString("vector", vstr)) {
      PmGeometry *geom;
      PmGeometryVector *vgeom;
      float dp;
      vector<PmVector3> vpts;
      pmSystem.getGeometry(vstr, &geom);

      if (!geom) {
        pm_ErrorReport (PM, "no vector named \"%s\" ", "*", vstr.c_str());
        return;
        }

      vgeom = dynamic_cast<PmGeometryVector*>(geom);
      vgeom->getCurrVector(vec);
      vgeom->getCurrPoints(vpts);

      if (dlist.getData("point", pdata1)) {
        if (!pm_CmdProcPosition(pdata1, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          return;
          }
        }
      else {
        pm_ErrorReport (PM, "no point given.", "*");
        return;
        }

      fprintf (stderr, ">>> point = %f %f %f \n", point[0], point[1], point[2]);
      points.push_back(point);
      vec.normalize();
      v = point - vpts[0];
      dp = v * vec;
      point = vpts[0] + dp*vec;
      points.push_back(point);
      v = points[1] - points[0];
      dist = v.length();

      pdata1.getString(pstr1);
      fprintf (stderr, ">>> distance between \"%s\" and \"%s\" = %f \n\n", vstr.c_str(),
              pstr1.c_str(), dist);
      }

    else if (dlist.getString("vector1", vstr)) {
      PmGeometry *geom1, *geom2, *geom3;
      PmGeometryVector *vgeom1, *vgeom2, *vgeom3;
      vector<PmVector3> vpts1, vpts2;
      PmVector3 pt1, pt2, v;
      int i1, i2;

      pmSystem.getGeometry(vstr, &geom1);
      
      if (!geom1) {
        pm_ErrorReport (PM, "no vector named \"%s\" ", "*", vstr.c_str());
        return;
        }

      vgeom1 = dynamic_cast<PmGeometryVector*>(geom1);
      vgeom1->getCurrPoints(vpts1);

      if (!dlist.getString("vector2", vstr)) {
        pm_ErrorReport (PM, "no vector2 given", "*");
        return;
        }

      pmSystem.getGeometry(vstr, &geom2);

      if (!geom2) {
        pm_ErrorReport (PM, "no vector named \"%s\" ", "*", vstr.c_str());
        return;
        }

      vgeom2 = dynamic_cast<PmGeometryVector*>(geom2);
      vgeom2->getCurrPoints(vpts2);

      if (!dlist.getInt("point1", i1)) {
        pm_ErrorReport (PM, "no point1 given", "*");
        return;
        }

      if ((i1 < 1) || (i1 > 2)) {
        pm_ErrorReport (PM, "point1 must be 1 or 2.", "*");
        return;
        }

      if (!dlist.getInt("point2", i2)) {
        pm_ErrorReport (PM, "no point2 given", "*");
        return;
        }

      if ((i2 < 1) || (i2 > 2)) {
        pm_ErrorReport (PM, "point2 must be 1 or 2.", "*");
        return;
        }

      if (dlist.getString("direction", vstr)) {
        pmSystem.getGeometry(vstr, &geom3);

        if (!geom3) {
          pm_ErrorReport (PM, "no vector named \"%s\" ", "*", vstr.c_str());
          return;
          }

        vgeom3 = dynamic_cast<PmGeometryVector*>(geom3);
        vgeom3->getCurrVector(dir);
        dir.normalize();
        in_dir = true;
        }

      pt1 = vpts1[i1-1];
      pt2 = vpts2[i2-1];
      v = pt1 - pt2;
      dist = v.length();
      //fprintf (stderr, ">>> pt1=%f %f %f\n", pt1[0], pt1[1], pt1[2]); 
      //fprintf (stderr, ">>> pt2=%f %f %f\n", pt2[0], pt2[1], pt2[2]); 

      if (in_dir) {
        float d = v*dir;
        fprintf (stderr, ">>> distance=%f  distance in direction=%f \n", dist, d); 
        }
      else {
        fprintf (stderr, ">>> distance=%f  value=%f %f %f\n", dist, v[0], v[1], v[2]); 
        }

      points.push_back(pt1);
      points.push_back(pt2);
      }

    // process picking two points //

    else if (dlist.getBoolean("pick", pick)) {
      fprintf (stderr, ">>> select two points \n"); 

      if (pick) {
        string cmd;
        stringstream dss;
        PmGraphicsInterface *grint = pmSystem.getGraphicsInterace();

        dss << "measure  distance ";
        cmd = dss.str();
        dss.str(std::string());
        grint->requestPickPoints(2, cmd);
        fprintf (stderr, ">>> select two points \n"); 
        }

      return;
      }

    else if (dlist.getData("point1", pdata1)) {
      if (!pm_CmdProcPosition(pdata1, point)) {
        pm_ErrorReport (PM, "couldn't process point data.", "*");
        return;
        }

      fprintf (stderr, ">>> point1=%f %f %f \n", point[0], point[1], point[2]);
      points.push_back(point);

      if (dlist.getData("point2", pdata2)) {
        if (!pm_CmdProcPosition(pdata2, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          return;
          }

        fprintf (stderr, ">>> point2=%f %f %f \n", point[0], point[1], point[2]);
        points.push_back(point);
        }

      if (points.size() != 2) {
        pm_ErrorReport (PM, "two points must be given.", "*");
        return;
        }

      v = points[1] - points[0];

      if (dlist.getVector("direction", dir)) {
        dist = dir * v;
        }
   
      else if (dlist.getString("vector", vstr)) {
        PmGeometry *geom;
        PmGeometryVector *vgeom;
        float dp;
        vector<PmVector3> vpts;
        pmSystem.getGeometry(vstr, &geom);

        if (!geom) {
          pm_ErrorReport (PM, "no vector named \"%s\" ", "*", vstr.c_str());
          return;
          }

        vgeom = dynamic_cast<PmGeometryVector*>(geom);
        vgeom->getCurrVector(vec);
        vec.normalize();
        dist = vec * v;
        }
      else {
        dist = v.length();
        }

      pdata1.getString(pstr1);
      pdata2.getString(pstr2);
      fprintf (stderr, ">>> distance=%f \n\n", dist); 
      }

    string geom_name;
    geom_name = "dist line";
    PmVector3 verts[2];
    verts[0] = points[0];
    verts[1] = points[1];

    if (!dist_geom) { 
      dist_geom = new PmGraphicsLine(geom_name, 2, verts);
      color.set(1,1,0);
      atts.setColor(color);
      atts.setDisjoint(true);
      atts.setLineWidth(2.0);
      dist_geom->setAttributes(atts);
      dist_geom->display();
      }
    else {
      dist_geom->update(2, verts);
      }

    if (fp) {
      string res;
      PmVector3 diff;
      float dist;

      diff = points[0] - points[1];

      if (in_dir) {
        dist = diff*dir;
        }
      else {
        dist = diff[2]; 
        }

      if (dlist.getString("residue", res)) {
        fprintf (fp, "%s %g\n", res.c_str(), dist);
        }
      else {
        fprintf (fp, "%g\n", dist);
        }

      fflush(fp);
      }
    }

  //===== rmsd between two domains =====//

  else if (data.name == "rmsd") {
    string seq1, seq2;
    string bname1, bname2;
    PmBody *body1, *body2;
    vector<PmVector3> coords1, coords2;

    //----- process body1 and sequence -----//

    if (!dlist.getString("body1", bname1)) {
      pm_ErrorReport (PM, "no body1 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    pmSystem.getBody (bname1, &body1);

    if (!body1) {
      pm_ErrorReport (PM, "no body names \"%s\".", "*", bname1.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (!dlist.getString("sequence1", seq1)) {
      pm_ErrorReport (PM, "no sequence1 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    body1->getDomainCoords(seq1, coords1);

    if (!coords1.size()) {
      pm_ErrorReport (PM, "could not process sequence1=%s.", "*", seq1.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    //----- process body2 and sequence -----//

    if (!dlist.getString("body2", bname2)) {
      pm_ErrorReport (PM, "no body2 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    pmSystem.getBody (bname2, &body2);

    if (!body2) {
      pm_ErrorReport (PM, "no body2 named \"%s\".", "*", bname2.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (!dlist.getString("sequence2", seq2)) {
      pm_ErrorReport (PM, "no sequence2 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    body2->getDomainCoords(seq2, coords2);

    if (!coords2.size()) {
      pm_ErrorReport (PM, "could not process sequence2=%s.", "*", seq2.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (coords1.size() != coords2.size()) {
      pm_ErrorReport (PM, "sequences are not equal: sequence=%s has %d atoms; sequnce2=%s has %d atoms.", "*", seq1.c_str(), coords1.size(), seq2.c_str(), coords2.size());
      pmSystem.setCmdError(true);
      return;
      }

    float err; 
    vector<float> wt;
    PmXform xform;

    pm_MathFitRms (coords1, coords2, wt, xform, err);
    fprintf (stderr, ">>> rmsd using %d coordinates \n", coords1.size());
    fprintf (stderr, ">>> rmsd=%g\n", err);
    }
  }

//*============================================================*
//*==========          pm_CmdMeasurement             ==========*
//*============================================================*
// process measurement command.

void
pm_CmdMeasurement (PmCmdDataList& dlist)
  {
  PmCmdData data, pdata;
  string dv, name, dname, type_str, bname, vname; 
  PmMeasurement *msr;
  PmMeasurementType mtype;

  vector<bool> global_frames;
  vector<PmVector3> points;
  bool bval;
  PmVector3 point;
  PmBody *body;
  vector<void*> pobjs;

  dlist.getNext (data);

  // create a new measurement and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmMeasurement::convMeasurementType(type_str, mtype);

      if (mtype == PM_MEASUREMENT_UNKNOWN) {
        pm_ErrorReport (PM, "unknown measurement type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no measurement type specified.", "*");
      return;
      }

    // create a measurement of a particular type and add it to the system //

    msr = PmMeasurement::create(name, mtype);
    pmSystem.addMeasurement(msr);

    if (mtype == PM_MEASUREMENT_DISTANCE) {
      PmDistMeasurement *dmsr = dynamic_cast<PmDistMeasurement*>(msr);

      if (dlist.getBoolean("global_frame1", bval)) {
        global_frames.push_back(bval);
        }

      if (dlist.getBoolean("global_frame2", bval)) {
        global_frames.push_back(bval);
        }

      if (dlist.getData("point1", pdata)) {
        if (!pm_CmdProcPosition(pdata, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          return;
          }
        //fprintf (stderr, ">>> point1 = %f %f %f \n", point[0], point[1], point[2]);
        points.push_back(point);
        }

      if (dlist.getString("body1", bname)) {
        pmSystem.getBody (bname, &body);

        if (!body) {
          pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
          return;
          }
        pobjs.push_back(body);
        }

      if (dlist.getData("point2", pdata)) {
        if (!pm_CmdProcPosition(pdata, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          return;
          }
        //fprintf (stderr, ">>> point2 = %f %f %f \n", point[0], point[1], point[2]);
        points.push_back(point);
        }

      if (dlist.getString("body2", bname)) {
        pmSystem.getBody (bname, &body);

        if (!body) {
          pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
          return;
          }

        pobjs.push_back(body);
        }

      if (global_frames.size() != 0) {
        dmsr->setGlobalFrames(global_frames);
        }

      if (points.size() != 2) {
        pm_ErrorReport (PM, "two points must be given.", "*");
        return;
        }
      else {
        dmsr->setPoints(points);
        }

      if (pobjs.size() != 2) {
        pm_ErrorReport (PM, "two bodies must be given.", "*");
        return;
        }
      else {
        dmsr->setObjs(pobjs);
        }
      }

    // process angle measurement //

    else if (mtype == PM_MEASUREMENT_ANGLE) {
      char bstr[20], pstr[20];
      PmAngleMeasurement *amsr = dynamic_cast<PmAngleMeasurement*>(msr);

      // angle between two vectors //

      if (dlist.getString("vector1", vname)) {
        PmGeometryVector *vgeom1, *vgeom2;
        PmGeometry *geom; 
        PmBody *body1, *body2;

        pmSystem.getGeometry(vname, &geom);
        vgeom1 = dynamic_cast<PmGeometryVector*>(geom);

        if (!geom || !vgeom1) {
          pm_ErrorReport (PM, "no vector named \"%s\".", "*", vname.c_str());
          return;
          }

        if (!dlist.getString("vector2", vname)) {
          pm_ErrorReport (PM, "no vector2 given.", "*");
          return;
          }
       
        pmSystem.getGeometry(vname, &geom);
        vgeom2 = dynamic_cast<PmGeometryVector*>(geom);

        if (!geom || !vgeom2) {
          pm_ErrorReport (PM, "no vector named \"%s\".", "*", vname.c_str());
          return;
          }

        if (dlist.getString("body1", bname)) {
          pmSystem.getBody (bname, &body1);

          if (!body1) {
            pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
            return;
            }
          }

        if (dlist.getString("body2", bname)) {
          pmSystem.getBody (bname, &body2);

          if (!body2) {
            pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
            return;
            }
          }

        amsr->setVectors(vgeom1, body1, vgeom2, body2);
        }

      // angle between two vectors givent by three points //

      else {
        for (int i = 1; i <= 3; i++) {
          sprintf (bstr, "body%d", i);
          sprintf (pstr, "point%d", i);

          if (dlist.getData(pstr, pdata)) {
            if (!pm_CmdProcPosition(pdata, point)) {
              pm_ErrorReport (PM, "couldn't process point data.", "*");
              return;
              }
            points.push_back(point);
            }
          else {
            pm_ErrorReport (PM, "no point \"%s\" given.", "*", pstr);
            return;
            }

          if (dlist.getString(bstr, bname)) {
            pmSystem.getBody (bname, &body);

            if (!body) {
              pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
              return;
              }
            pobjs.push_back(body);
            }
          else {
            pm_ErrorReport (PM, "no body \"%s\" given.", "*", bstr);
            return;
            }
          }

        if (points.size() != 3) {
          pm_ErrorReport (PM, "three points must be given.", "*");
          return;
          }
        else {
          amsr->setPoints(points);
          }

        if (pobjs.size() != 3) {
          pm_ErrorReport (PM, "three bodies must be given.", "*");
          return;
          }
        else {
          amsr->setObjs(pobjs);
          }
        }
      }

    // process position measurement //

    else if (mtype == PM_MEASUREMENT_POSITION) {
      PmPositionMeasurement *pmsr = dynamic_cast<PmPositionMeasurement*>(msr);

      if (dlist.getData("point", pdata)) {
        if (!pm_CmdProcPosition(pdata, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          return;
          }
        points.push_back(point);
        }

      if (dlist.getString("body", bname)) {
        pmSystem.getBody (bname, &body);

        if (!body) {
          pm_ErrorReport (PM, "no body named \"%s\". ", "*", bname.c_str());
          return;
          }
        pobjs.push_back(body);
        }

      if (points.size() != 1) {
        pm_ErrorReport (PM, "a point must be given.", "*");
        return;
        }
      else {
        pmsr->setPoints(points);
        }

      if (pobjs.size() != 1) {
        pm_ErrorReport (PM, "a body must be given.", "*");
        return;
        }
      else {
        pmsr->setObjs(pobjs);
        }
      }
    }

  else {
    pmSystem.getMeasurement(data.name, &msr);
    }

  if (!msr) {
    pm_ErrorReport (PM, "no measurement named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process measurement parameters //

  PmVector3 color;
  bool show_set = true;

  while (dlist.getNext(data)) {
    if (data.name == "color") {
      data.getVector(color);
      msr->setColor(color);
      }
    }

  if (show_set) {
    msr->display ();
    }
  }

