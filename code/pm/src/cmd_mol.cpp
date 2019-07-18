
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
//*==========              pm_CmdMoleculeRead        ==========*
//*============================================================*
// process molecule read command.

void
pm_CmdMoleculeRead (PmMolecule *mol, PmCmdDataList& dlist)
  {
  string dv, format_str, fname;
  vector<PmMode> modes;
  PmDbInterfaceSelect db_sel;
  PmDbType format;
  PmCmdData data;

  dlist.getNext (data);

  if (data.name == "modes") {
    if (!dlist.getString("format", format_str)) {
      pm_ErrorReport (PM, "no format specified.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    if (!dlist.getString("file", fname)) {
      pm_ErrorReport (PM, "no file name specified.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    format = pm_DbFormatConv(format_str);

    if (format == PM_DB_UNKNOWN) {
      pm_ErrorReport (PM, "unknown format type \"%s\".", "*", format_str.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    PmDbInterface *db = db_sel.create("modes", format);

    if (!db) {
      pm_ErrorReport ("pm> ", "unknown format \"%s\".", "*", format_str.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    db->open(fname, PM_DB_MODE_READ, "");
    db->getModes(modes);

    if (modes.size()) {
      mol->setModes(modes);
      }
    }
  }

//*============================================================*
//*==========              pm_CmdMolecule            ==========*
//*============================================================*
// process molecule command.                                 

void
pm_CmdMolecule (PmCmdDataList& dlist)
  {
  int model;
  string dv, dbname, name, type, model_name; 
  PmMolecule *mol;
  PmMoleculeType mol_type;
  PmExtent extent;
  PmDbInterface *db;
  PmCmdData data;

  dlist.getNext (data);

  if (data.name == "read") {
    if (!dlist.getString("name", name)) {
      pm_ErrorReport (PM, "no molecule name given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    if (!dlist.getString ("database", dbname)) {
      pm_ErrorReport (PM, "no database given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    dlist.getString ("type", type);

    if (type != "") { 
      PmMolecule::convMoleculeType(type, mol_type);

      if (type == "rna") {
        mol_type = PM_MOLECULE_NUCLEIC_ACID; 
        }
      else if (type == "dna") {
        mol_type = PM_MOLECULE_NUCLEIC_ACID; 
        }
      else { 
        PmMolecule::convMoleculeType(type, mol_type);
        }

      if (mol_type == PM_MOLECULE_UNKNOWN) {
        pm_ErrorReport (PM, "unknown molecule type \"%s\".", "*", type.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else {
      type = "protein";
      mol_type = PM_MOLECULE_PROTEIN; 
      }

    pmSystem.getDatabase(dbname, &db);

    if (!db) {
      pm_ErrorReport (PM, "no database named \"%s\".", "*", dbname.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (dlist.getString("model_name", model_name)) {
      mol = db->getMolecule(model_name, mol_type);

      if (!mol) {
        pm_ErrorReport (PM, "could not read molecule for type \"%s\" and model name %s.", 
                        "*", type.c_str(), model_name.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }

    else if (dlist.getInt("model", model)) {
      mol = db->getMolecule(model, mol_type);

      if (!mol) {
        pm_ErrorReport (PM, "could not read molecule for type \"%s\" and model %d.", 
                        "*", type.c_str(), model);
        pmSystem.setCmdError(true);
        return;
        }
      }

    else {
      mol = db->getMolecule(1, mol_type);

      if (!mol) {
        pm_ErrorReport (PM, "could not read molecule for type \"%s\" and model %d.", 
                        "*", type.c_str(), 1);
        pmSystem.setCmdError(true);
        return;
        }
      }

    pm_PrintMsg (CMDPR, "number of molecule atoms = [%d]", mol->getNumAtoms());
    mol->setName(name);
    mol->setType(mol_type);
    pm_PrintMsg (CMDPR, "molecule type = %s ", type.c_str()); 

    // add molecule to system //
    pmSystem.addMolecule (mol);
    mol->getExtent(extent);
   
    pm_PrintMsg (CMDPR, "extent  min (%g %g %g)  max (%g %g %g) ", 
             extent.min[0], extent.min[1], extent.min[2],
             extent.max[0], extent.max[1], extent.max[2]);

    
    PmVector3 center; 
    center = (extent.min + extent.max) / 2.0;
    pm_PrintMsg (CMDPR, "center: (%g %g %g)  ", center[0], center[1], center[2] );

    // set extent //
    //pmSystem.setExtent (extent);
    return;
    }

  else if (data.name == "set") {
    dlist.getString("name", name);
    pmSystem.getMolecule(name, &mol);
    pmSystem.setCurrentMolecule(mol);
    return;
    }

  pmSystem.getMolecule (data.name, &mol);

  if (!mol) {
    pm_ErrorReport (PM, "no molecule named \"%s\".", "*", data.name.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  // a lot of commands for molecules and domains are the same //
  // so call the domain command processing functions.         //

  while (dlist.getNext(data, dv)) {

    // ===== atoms ===== //

    if (data.name == "atoms") {
      pm_CmdDomainAtoms (mol, dlist);
      }

    // ===== backbone ===== //

    else if (data.name == "backbone") {
      pm_CmdDomainBackbone (mol, dlist);
      }

    // ===== bonds ===== //

    else if (data.name == "bonds") {
      pm_CmdDomainBonds (mol, dlist);
      }

    // ===== fit ===== //

    else if (data.name == "fit") {
      pm_CmdDomainBonds (mol, dlist);
      }

    // ===== print ===== //

    else if (data.name == "print") {
      mol->print(); 

      if (dv == "secondary_structure") {
        mol->printSecondaryStructure(); 
        }
      }

    // ===== read ===== //

    else if (data.name == "read") {
      pm_CmdMoleculeRead(mol, dlist);
      return;
      }

    // ===== surface ===== //

    else if (data.name == "surface") {
      pm_CmdDomainSurface (mol, dlist);
      }

    // ===== xform ===== //

    else if (data.name == "xform") {
      pm_CmdDomainXform (mol, dlist);
      }
    }
  }

