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

///////////////////////////////////////////////////////////////
//               P D B    I n t e r f a c e                 //
/////////////////////////////////////////////////////////////

// define pdb file processing table

typedef void (PmDbPdbInterface::*Pdbfp)(char*);

typedef struct PmDbPdbProc {
  char *name;
  Pdbfp func; 
  } PmDbPdbProc;

static PmDbPdbProc  pdb_table[] = {
    {"atom",     &PmDbPdbInterface::pm_DbPdbAtomProc},
    {"conect",   &PmDbPdbInterface::pm_DbPdbConectProc},
    {"endmdl",   &PmDbPdbInterface::pm_DbPdbEndModelProc},
    {"helix",    &PmDbPdbInterface::pm_DbPdbHelixProc},
    {"hetatm",   &PmDbPdbInterface::pm_DbPdbHetAtomProc},
    {"model",    &PmDbPdbInterface::pm_DbPdbModelProc},
    {"seqres",   &PmDbPdbInterface::pm_DbPdbSeqresProc},
    {"sheet",    &PmDbPdbInterface::pm_DbPdbSheetProc},
    {NULL,       NULL}};

static char
*atom_fmt = "%6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f           %2s \n";

static char
*model_fmt = "MODEL     %4d   %s\n";


//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmDbPdbInterface::PmDbPdbInterface(const string name) { 
  this->name = name;
  molecule_type = PM_MOLECULE_UNKNOWN;
  protein = NULL;
  water = NULL;
  nucleotide = NULL;
  misc = NULL;
  use_misc = false;
  model_id = 1;
  }

void
PmDbPdbInterface::getCurrentMolecule(PmMoleculeType type, PmMolecule **mol)
  {
  //fprintf (stderr, "    >>> get curr mol  type [%d] \n", type); 
  char str[200];
  memset (str, 0, 200);
  sprintf (str, "%s:%d", molecule_name.c_str(), model_id);
  string mol_name;
  mol_name = str;
  PmMoleculeType mtype;

  if (molecule_type != PM_MOLECULE_UNKNOWN) {
    mtype = molecule_type; 
    }
  else {
    mtype = type; 
    }

  if ((mtype == PM_MOLECULE_PROTEIN) || (mtype == PM_MOLECULE_AMINO_ACID)) {
    if (!protein) {
      protein = new PmMolecule(mol_name);
      protein->setModelName(model_name);
      protein_molecules.push_back(protein);
      model_names.push_back(model_name);
      }

    *mol = protein;
    }

  else if (mtype == PM_MOLECULE_NUCLEOTIDE) { 
    if (!nucleotide) {
      nucleotide = new PmMolecule(mol_name);
      nucleotide_molecules.push_back(nucleotide);
      }

    *mol = nucleotide;
    }

  else if (mtype == PM_MOLECULE_WATER) { 
    if (!water) {
      water = new PmMolecule(mol_name);
      water_molecules.push_back(water);
      }

    *mol = water;
    }

  else { 
    if (!misc) {
      misc = new PmMolecule(mol_name);
      misc_molecules.push_back(misc);
      }

    *mol = misc;
    }
  }

//*============================================================*
//*==========              open                      ==========*
//*============================================================*
// open pdb file.

void
PmDbPdbInterface::open (const string name, PmDbModeType otype, string dtype) 
  {
  bool verbose = pmSystem.getCmdVerbose();

  if (verbose) {
    fprintf (stderr, "\n\n>>>>>> PDB: open file=%s", name.c_str());
    }

  if (otype == PM_DB_MODE_READ) {
    if (name.find(".gz") != string::npos) {
      string cmd;
      cmd = "zcat -d " + name;
      fp = popen (cmd.c_str(), "r");
      //fprintf (stderr, ">>> open compressed file. \n");
      }
    else {
      fp = fopen (name.c_str(), "r");
      }

    if (!fp) {
      pm_ErrorReport (PM, "can't open \"%s\" for reading.", "*", name.c_str());
      return;
      }

    if (verbose) {
      fprintf (stderr, " for reading. \n");
      }

    if (dtype != "") {
      PmMolecule::convMoleculeType(dtype, molecule_type);
      fprintf (stderr, ">>> data is %s \n", dtype.c_str());
      }

    this->read(name);
    this->close();
    }
  else if (otype == PM_DB_MODE_APPEND) {
    fp = fopen (name.c_str(), "a");

    if (!fp) {
      pm_ErrorReport (PM, "can't open \"%s\" for appending.", "*", name.c_str());
      return;
      }

    if (verbose) {
      fprintf (stderr, " for appending. \n");
      }
    }
  else {
    fp = fopen (name.c_str(), "w");

    if (!fp) {
      pm_ErrorReport (PM, "can't open \"%s\" for writing.", "*", name.c_str());
      return;
      }

    if (verbose) {
      fprintf (stderr, " for writing. \n");
      }
    }
  }

//*============================================================*
//*==========              close                     ==========*
//*============================================================*

void
PmDbPdbInterface::close() {
  if (fp) {
    fclose (fp);
    }
  }

//*============================================================*
//*==========              getMolecule               ==========*
//*============================================================*
// get molecule.      

PmMolecule *
PmDbPdbInterface::getMolecule(int model, PmMoleculeType type) 
  {
  //fprintf (stderr, ">>>>>> PmDbPdbInterface::getMolecule \n");
  //fprintf (stderr, ">>> model=%d  type=%d \n", model, type);
  PmMolecule *mol = NULL;

  if ((model < 1) || (model > model_id)) {
    return mol;
    }

  if ((type == PM_MOLECULE_PROTEIN) && protein_molecules.size()) { 
    mol = protein_molecules[model-1];
    }

  else if ((type == PM_MOLECULE_NUCLEOTIDE) && nucleotide_molecules.size()) { 
    mol = nucleotide_molecules[model-1];
    }

  else if ((type == PM_MOLECULE_NUCLEIC_ACID) && nucleotide_molecules.size()) { 
    mol = nucleotide_molecules[model-1];
    }

  else if (misc_molecules.size()) { 
    mol = misc_molecules[model-1];
    }

  return mol;
  }

//===== get molecule by name =====//

PmMolecule *
PmDbPdbInterface::getMolecule(string model, PmMoleculeType type) 
  {
  //fprintf (stderr, ">>>>>> PmDbPdbInterface::getMolecule \n");
  //fprintf (stderr, ">>> model=%s  type=%d \n", model.c_str(), type);
  PmMolecule *mol = NULL;
  string model_name;

  if ((type == PM_MOLECULE_PROTEIN) && protein_molecules.size()) { 
    for (unsigned int i = 0; i < protein_molecules.size(); i++) {
      protein_molecules[i]->getModelName(model_name);
      //fprintf (stderr, ">>> mol=%d  model_name=%s \n", i, model_name.c_str());

      if (model_name == model) {
        mol = protein_molecules[i];
        }
      }
    }

  return mol;
  }

//*============================================================*
//*==========              getModelNames             ==========*
//*============================================================*

void 
PmDbPdbInterface::getModelNames(vector<string>& model_names) 
  {
  model_names = this->model_names;
  }

//*============================================================*
//*==========              read                      ==========*
//*============================================================*
// read pdb file.

void
PmDbPdbInterface::read(const string mol_name) 
  {

  int num;
  char line[DB_MAX_LINE];

  if (!fp) {
    return;
    }

  bool verbose = pmSystem.getCmdVerbose();

  // parse the file //

  if (verbose) fprintf (stderr, "        read pdb file ... ");
  num = 0;
  molecule_name = mol_name;

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    //fprintf (stderr, "%d[%s]\n", num, line);
    num += 1;

    for (int i = 0; pdb_table[i].name != NULL; i++) {
      int n = strlen (pdb_table[i].name);

      if (!strncasecmp(line, pdb_table[i].name, n)) {
        (this->*pdb_table[i].func) (line);
        }
      }
    }

  if (verbose) {
    fprintf (stderr, "done.\n");
    fprintf (stderr, "\n");
    fprintf (stderr, "        read %d protein model(s)\n", protein_molecules.size());
    fprintf (stderr, "             %d nucleotide model(s)\n", nucleotide_molecules.size());
    fprintf (stderr, "             %d misc model(s) \n", misc_molecules.size());
    fprintf (stderr, "\n");
    }

  for (int j = 0; j < 2; j++) {
    PmMolecule *mol = NULL; 

    if ((j == 0) && protein_molecules.size()) {
      mol = protein_molecules[0];
      if (verbose) fprintf (stderr, "\n    ============= protein ============\n");
      }

    else if ((j == 1) && nucleotide_molecules.size()) {
      mol = nucleotide_molecules[0];
      if (verbose) fprintf (stderr, "\n    ============ nucleotide ============\n");
      }

    if (mol) {
      int num_atoms = mol->getNumAtoms();
      int num_chains = mol->getNumChains();
      if (verbose) fprintf (stderr, "    >>> %d atoms   %d chains \n", num_atoms, num_chains);

      vector<PmChain*> clist;
      mol->getChains(clist);

      if (verbose) {
        fprintf (stderr, "    --------- chains --------- \n");

        for (unsigned int i = 0; i < clist.size(); i++) {
          fprintf (stderr, "    %d  id = %c  number of residues = %d\n", i+1, clist[i]->id,
                   clist[i]->num_residues);
          }
        }
      }
    }
  }

//*============================================================*
//*==========              writeHeader               ==========*
//*============================================================*
// write a header.

void
PmDbPdbInterface::writeHeader()
  {
  fprintf (this->fp, "REMARK  Protein Mechanica simulation domain file \n");
  fflush (this->fp);
  }

//*============================================================*
//*==========              pm_DbPdbOutputAtom        ==========*
//*============================================================*
// write an atom. 

static void
pm_DbPdbOutputAtom (FILE *fp, PmAtom *atom, int atom_num, char chain_id, int seq)
  {

  char elem_name[3], cmpd_name[3];
  PmAtomName atom_name, aname;
  int n, m;

  strcpy (atom_name, atom->name);
  PmAtom::convType(atom->type, elem_name);
  pm_CmpdGetName(atom->ctype, cmpd_name);

  float x, y, z, mass = 1.0;
  x = 10.0 * atom->pos[0];
  y = 10.0 * atom->pos[1];
  z = 10.0 * atom->pos[2];

  for (int k = 0; k < 3; k++) {
    cmpd_name[k] = toupper(cmpd_name[k]);
    }

  // pad atom name with blanks //

  n = sizeof(PmAtomName);
  m = strlen(atom_name);

  for (int k = 0; k < n-1; k++) {
    aname[k] = ' ';
    }

  aname[n-1] = '\0';

  for (int k = 0; k < m; k++) {
    aname[k+1] = atom_name[k];
    }

  //fprintf (stderr, "12345   n = %d   m = %d \n", n, m);
  //fprintf (stderr, "%s\n", aname);

  fprintf (fp, atom_fmt, "ATOM  ", atom_num, aname, ' ', cmpd_name, chain_id,
           seq, ' ', x, y, z, 1.0, mass, elem_name);
  //fflush (fp);`
  }

//*============================================================*
//*==========              writeDynDomainAtoms       ==========*
//*============================================================*
// write the atoms for a list of dynamic domains.

void
PmDbPdbInterface::writeDynDomainAtoms(const int step, const bool init,
                                      vector<PmMolecule*> domains, 
                                      vector<PmXform> xforms)
  {
  #define ndbg_PmDbPdbInterface_writeDynDomainAtoms

  #ifdef dbg_PmDbPdbInterface_writeDynDomainAtoms
  fprintf (stderr, ">>>>>> PmDbPdbInterface::writeDynDomainAtoms step = %d\n", 
           step);
  fprintf (stderr, "   >>> num domains = %d\n", domains.size());
  #endif

  vector<PmAtom> atoms;
  PmMatrix3x3 mat;
  PmVector3 pos, xpos, center, translation;
  PmXform xform;
  vector<char> chains; 
  vector<PmChain*> clist;
  char chain_id;
  bool chain_found;
  PmAtomFilter filter;
  string dname;

  // get list of chains //

  for (unsigned int i = 0; i < domains.size(); i++) {
    domains[i]->getChains(clist);
    #ifdef dbg_PmDbPdbInterface_writeDynDomainAtoms
    domains[i]->getName(dname);
    //fprintf (stderr, "   >>> domain \"%s\" \n", dname.c_str());
    //fprintf (stderr, "       clist.size = %d \n", clist.size());
    #endif 

    for (unsigned int j = 0; j < clist.size(); j++) {
      chain_id = clist[j]->id;
      chain_found = false;

      for (unsigned int k = 0; k < chains.size(); k++) {
        if (chain_id == chains[k]) {
          chain_found = true;
          }
        }

      if (!chain_found) {
        chains.push_back(chain_id);
        //fprintf (stderr, "       push chain = %c \n", chain_id);
        }
      }
    }


  // write atoms //

  #ifdef dbg_PmDbPdbInterface_writeDynDomainAtoms
  //fprintf (stderr, "--- num domains=%d --- \n", domains.size()); 
  //for (unsigned int i = 0; i < domains.size(); i++) {
    //string dname;
    //domains[i]->getName(dname);
    //fprintf (stderr, "domain %s \n", dname.c_str()); 
    //}

  //fprintf (stderr, "---- write atoms ------- \n");
  #endif

  fprintf (fp, "REMARK  step = %d \n", step);

  for (unsigned int k = 0; k < chains.size(); k++) {
    chain_id = chains[k]; 

    for (unsigned int i = 0; i < domains.size(); i++) {
      domains[i]->getChainAtoms(chain_id, atoms);
      string dname;
      domains[i]->getName(dname);
      //fprintf (stderr, "domain %s \n", dname.c_str()); 

      if (atoms.size()) {
        #ifdef dbg_PmDbPdbInterface_writeDynDomainAtoms
        fprintf (stderr, ">>> domain=%s chain=%c  seq=%d-%d\n", dname.c_str(),
                          chain_id, atoms[0].seq, atoms[atoms.size()-1].seq);
        #endif
        mat = xforms[i].matrix;
        center = xforms[i].center;
        translation = xforms[i].translation;
        //domains[i]->getAtoms(atoms);

        for (unsigned int j = 0; j < atoms.size(); j++) {
          pos = atoms[j].pos;
          xpos = mat*(pos-center) + translation + center;
          atoms[j].pos = xpos;
          pm_DbPdbOutputAtom (fp, &atoms[j], atoms[j].id, atoms[j].chain, atoms[j].seq);
          }
        }
      }
    }

  fprintf (fp, "END\n");
  fflush(fp);
  }

//*============================================================*
//*==========              writeDomainAtoms          ==========*
//*============================================================*
// write the atoms for a list of domains.

void 
PmDbPdbInterface::writeDomainAtoms(vector<PmMolecule*> domains, vector<string> chains,
                                   string one_chain, bool no_hydrogen, string start_res, 
                                   string res_names, int model, string model_name)
  {
  //fprintf (stderr, "\n>>>>>> PmDbPdbInterface::writeDomainAtoms \n");

  if (!fp) {
    return;
    }

  int num_dom = domains.size();
  //fprintf (stderr, ">>> num domains=%d. \n", num_dom);
  //fprintf (stderr, ">>> no hydrogen=%d. \n", no_hydrogen);
  vector<PmAtom> atoms;
  float mass;
  int start_res_num, res_count;
  char output_chain_id, chain_id, start_res_chain;
  mass = 1.0;
  res_count = 0;
  //fprintf (stderr, ">>> chains=%d. \n", chains.size());
  //fprintf (stderr, ">>> one chain=\"%s\" \n", one_chain.c_str());

  if (start_res.size()) {
    sscanf (start_res.c_str(), "%c%d", &start_res_chain, &start_res_num);
    //fprintf (stderr, ">>> start res chain=%c\n", start_res_chain); 
    //fprintf (stderr, ">>> start res num=%d\n", start_res_num); 
    //fprintf (stderr, ">>> num res names=%d\n", res_names.size()); 
    }
  else {
    start_res_chain = 0;
    start_res_num = 0;
    }

  if (model > 0) {
    fprintf (fp, model_fmt, model, model_name.c_str());
    }

  for (int i = 0; i < num_dom; i++) {
    PmMolecule *dom = domains[i];
    dom->getAtoms(atoms);
    int num_atoms = atoms.size();
    //fprintf (stderr, ">>> num atoms=%d \n", num_atoms); 
    int atom_num = 1;
    int seq;

    // output atoms as a single chain  //

    if (one_chain != "") {
      int useq = -1;
      bool use_seq = true;
      int start_seq;
      sscanf (one_chain.c_str(), "%c[%d]", &chain_id, &useq);
      //fprintf (stderr, ">>> chain[%c] seq[%d]. \n", chain_id, useq); 

      if (useq != -1) {
        use_seq = false;
        start_seq = atoms[0].seq;
        }

      for (int j = 0; j < num_atoms; j++) {
        if (use_seq) {
          seq = atoms[j].seq;
          }
        else {
          seq = useq + atoms[j].seq - start_seq;
          }

        if (no_hydrogen) {
          if (atoms[j].type != PM_ATOM_ELEMENT_HYDROGEN) {  
            pm_DbPdbOutputAtom(fp, &atoms[j], atom_num, chain_id, seq);
            atom_num += 1;
            }
          }
        else {
          pm_DbPdbOutputAtom(fp, &atoms[j], atom_num, chain_id, seq);
          atom_num += 1;
          }
        }
      }

    //===== output atoms with their original chain  =====//

    else if (!chains.size()) {
      //fprintf (stderr, ">>> start_res_chain=%c \n", start_res_chain); 
      atom_num = 1;
      int last_seq = -1;
      PmCompoundType ctype=PM_COMPOUND_UNKNOWN; 
      PmMoleculeType mtype;
      char cstr[2] = { '\0' }, cmpd_name[3];

      for (int j = 0; j < num_atoms; j++) {
        chain_id = atoms[j].chain;
        seq = atoms[j].seq;

        if ((chain_id == start_res_chain) && (last_seq != seq)) {
          last_seq = seq;

          if ((res_count < res_names.size()) && (seq == start_res_num+res_count)) {
            cstr[0] = tolower(res_names[res_count]);
            pm_CmpdGetName(atoms[j].ctype, cmpd_name);
            pm_CmpdConvType (cstr, ctype, mtype);
            fprintf (stderr, ">>> replace seq=%d  %s -> %s \n", seq, cmpd_name, cstr);
            res_count += 1;
            }
          }

        if (ctype != PM_COMPOUND_UNKNOWN) {
          atoms[j].ctype = ctype;
          }

        if (no_hydrogen) {
          if (atoms[j].type != PM_ATOM_ELEMENT_HYDROGEN) {
            pm_DbPdbOutputAtom(fp, &atoms[j], atom_num, chain_id, seq);
            atom_num += 1;
            }
          }
        else {
          pm_DbPdbOutputAtom(fp, &atoms[j], atom_num, chain_id, seq);
          atom_num += 1;
          }
        }
      }

    // output atoms with a given chain that may //
    // be modfied.                              //

    else {

      // get the residues to output //

      for (unsigned int k = 0; k < chains.size(); k++) {
        int start_seq; 
        int output_seq = -1;
        string src, dst; 
        bool use_seq = true;

        string s = chains[k];
        int n = s.length();
        size_t found = s.find('=');

        if (found != string::npos) {
          src.assign(s, 0, found); 
          dst.assign(s, found+1, n-found+1); 
          fprintf (stderr, ">>> source=%s  destination=%s\n", src.c_str(), 
                   dst.c_str()); 
          sscanf (dst.c_str(), "%c[%d]", &output_chain_id, &output_seq);
          fprintf (stderr, "    output chain=%c  output seq=%d\n", output_chain_id, 
                   output_seq);
          }
        else {
          src = s;
          chain_id = s[0];
          output_chain_id = chain_id;
          fprintf (stderr, ">>> source=%s\n", src.c_str());
          }

        vector<PmResidue*> rlist;
        PmResidue *res;
        dom->getResidues(src, rlist);
        fprintf (stderr, ">>> num residues=%d \n", rlist.size()); 

        if (!rlist.size()) {
          continue;
          }

        if (output_seq != -1) {
          use_seq = false;
          start_seq = rlist[0]->atoms[0]->seq;
          }

        for (unsigned int j = 0; j < rlist.size(); j++) {
          res = rlist[j];

          if ((output_chain_id == start_res_chain) && (res_count < res_names.size())) {
            if (res->atoms[k]->seq == start_res_num+res_count) {
              fprintf (stderr, ">>> replace seq=%d \n", res->atoms[k]->seq);
              res_count += 1;
              }
            }

          for (int k = 0; k < res->num_atoms; k++) {
            PmAtom *atom = res->atoms[k];

            if (use_seq) {
              seq = atom->seq;
              }
            else { 
              seq = output_seq + atom->seq - start_seq;
              }

            if (no_hydrogen) {
              if (atom->type != PM_ATOM_ELEMENT_HYDROGEN) {
                pm_DbPdbOutputAtom(fp, atom, atom_num, output_chain_id, seq);
                atom_num += 1;
                }
              }
            else {
              pm_DbPdbOutputAtom(fp, atom, atom_num, output_chain_id, seq);
              atom_num += 1;
              }
            }
          }
        }
      }
    }

  if (model > 0) {
    fprintf (fp, "ENDMDL\n");
    }

  fflush (fp);
  }

//*============================================================*
//*==========              writeDomainSequence       ==========*
//*============================================================*
// write the sequences for a list of domains.

void
PmDbPdbInterface::writeDomainSequence(vector<PmMolecule*> domains)
  {
  if (!fp) {
    return;
    }

  int n;
  int num_dom = domains.size();
  fprintf (stderr, "    >>> num domains[%d]. \n", num_dom);
  vector<PmResidue*> rlist;
  PmResidue *res;
  char cname[3];
  n = 0;

  for (int i = 0; i < num_dom; i++) {
    PmMolecule *dom = domains[i];
    dom->getResidues("", rlist);

    for (int j = 0; j < rlist.size(); j++) {
      res = rlist[j];
      pm_CmpdGetName(res->ctype, cname, 1);
      fprintf (fp, "%s", cname);
      n += 1;

      if (n == 40) {
        fprintf (fp, "\n");
        n = 0;
        }
      }
    }

  fprintf (fp, "\n");
  }

//*============================================================*
//*==========              setAtomNumber             ==========*
//*============================================================*

void 
PmDbPdbInterface::setAtomNumber(const char chain_id, const int num)
  {
  for (unsigned int i = 0; i < atom_numbers.size(); i++) {
    if (atom_numbers[i].chain_id == chain_id) {
      atom_numbers[i].start = num;
      return;
      }
    } 

  PmDbAtomSeq seq;
  seq.chain_id = chain_id;
  seq.start = num;
  atom_numbers.push_back(seq);
  }

void 
PmDbPdbInterface::getAtomNumber(const char chain_id, int& num)
  {
  for (unsigned int i = 0; i < atom_numbers.size(); i++) {
    if (atom_numbers[i].chain_id == chain_id) {
      num = atom_numbers[i].start;
      return;
      }
    }
 
  num = 1;
  }

//////////////////////////////////////////////////////////////////
//         functions to process pdb file entries               //
////////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_DbPdbAtomProc          ==========*
//*============================================================*
// process ATOM entry.

void
PmDbPdbInterface::pm_DbPdbAtomProc (char *line) 
  {

  int id, res_seq, n, i;
  char cname[4], str[10], chain_id;
  PmAtomName name;
  float occupancy, temp;
  PmVector3 pos;
  PmAtomElementType type;
  PmCompoundType ctype;
  PmMoleculeType mtype;
  PmMolecule *mol;

  //fprintf (stderr, "\n ----------  pm_DbPdbAtomProc  ---------- \n");
  //fprintf (stderr, " line [%s] \n", line);
  pm_DbRecordIntGet (line, 7, 11, id);
  pm_DbRecordStrGet (line, 13, 16, str);

  for (i = 0, n = 0; str[i]; i++) {
    if (str[i] != ' ') name[n++] = str[i];
    }
  name[n++] = '\0';

  pm_DbRecordStrGet (line, 18, 20, cname);
  //fprintf (stderr, " cname [%s] \n", cname);


  // alternate location indicator. i guess use the 1st one //

  if ((line[16] != 'A') && (line[16] != ' ')) {
    return;
    }

  pm_DbRecordIntGet (line, 23, 26, res_seq);
  pm_DbRecordRealGet (line, 31, 38, pos[0]);
  pm_DbRecordRealGet (line, 39, 46, pos[1]);
  pm_DbRecordRealGet (line, 47, 54, pos[2]);

  pm_DbRecordRealGet (line, 55, 60, occupancy);
  pm_DbRecordRealGet (line, 61, 66, temp);
  pm_DbRecordStrGet (line, 77, 78, str);

  if (isdigit(name[0])) {
    type = PmAtom::convType(&name[1]);
    }
  else {
    type = PmAtom::convType(&name[0]);
    }

  chain_id = line[21];

  if (chain_id == ' ') {
    chain_id = '1';
    } 

  if (res_seq == 0) {
    res_seq = 1;
    }

  // convert angstroms into nm //

  pos[0] /= 10.0;
  pos[1] /= 10.0;
  pos[2] /= 10.0;

  // determine compound and molecule type. use the  //
  // misc symbol to make sure atoms with a          //
  // a non-standard notation are not misclassified. //

  if (molecule_type == PM_MOLECULE_MISCELLANEOUS) {
    ctype = PM_COMPOUND_UNKNOWN;
    mtype = molecule_type;
    }
  else {
    pm_CmpdConvType (cname, ctype, mtype);
    }

  // add atom to molecule obj //

  getCurrentMolecule(mtype, &mol);
  mol->addAtom(id, name, type, pos, chain_id, ctype, mtype, res_seq);
  }

//*============================================================*
//*==========              pm_DbPdbConectProc        ==========*
//*============================================================*
// process ATOM entry.

void
PmDbPdbInterface::pm_DbPdbConectProc (char *line)
  {
  char *s;
  int n, conn[1000];

  s = strtok (line, " ");
  s = strtok (NULL, " ");
  n = 0;

  // store connect using original atom ids //

  while (s && (*s != '\n')) {
    conn[n++] = atoi(s);
    s = strtok (NULL, " ");
    }

  // NOTE: this needs to be fixed to check which molecule
  //       to store the connectivity.

  if (protein) { 
    protein->addConnect (n, conn);
    }

  else if (nucleotide) { 
    nucleotide->addConnect (n, conn);
    }

  else if (misc) { 
    misc->addConnect (n, conn);
    }
  }

//*============================================================*
//*==========              pm_DbPdbEndModelProc      ==========*
//*============================================================*
// process ENDMDL entry.

void
PmDbPdbInterface::pm_DbPdbEndModelProc(char *line) {
  protein = NULL;
  water = NULL;
  nucleotide = NULL;
  misc = NULL;
  }

//*============================================================*
//*==========          pm_DbPdbHetAtomProc           ==========*
//*============================================================*
// process HETATOM entry

void
PmDbPdbInterface::pm_DbPdbHetAtomProc (char *line)
  {
  char str[80];

  int id, res_seq;
  
  float x, y, z, temp;

  pm_DbRecordIntGet (line, 7, 11, id);
  pm_DbRecordStrGet (line, 13, 16, str);
  //atom.name = (char *)strdup(str);

  pm_DbRecordIntGet (line, 23, 26, res_seq);
  pm_DbRecordRealGet (line, 31, 38, x);
  pm_DbRecordRealGet (line, 39, 46, y);
  pm_DbRecordRealGet (line, 47, 54, z);

  pm_DbRecordRealGet (line, 61, 66, temp);
  pm_DbRecordStrGet (line, 77, 78, str);

  pm_DbPdbAtomProc (line);
  }

//*============================================================*
//*==========              pm_DbPdbModelProc         ==========*
//*============================================================*
// process MODEL entry.

void
PmDbPdbInterface::pm_DbPdbModelProc(char *line)
  {
  int id;
  char str[80];
  pm_DbRecordIntGet (line, 11, 14, id);
  model_id = id;
  pm_DbRecordStrGet (line, 18, 80, str);
  //fprintf (stderr, ">>> model id=%d \n", id);
  //fprintf (stderr, ">>> model name=%s \n", str);
  model_name.clear(); 

  for (int i = 0; str[i] != '\0'; i++) {
    if (str[i] != '\n') { 
      model_name.push_back(str[i]); 
      }
    }
  }

//*============================================================*
//*==========           pm_DbPdbSeqresProc           ==========*
//*============================================================*
// process HETATOM entry

void
PmDbPdbInterface::pm_DbPdbSeqresProc(char *line)
  {
  return;

  int res_seq, num_res, n;

  char chain_id, *s, str[80];

  //PmResName *residues;

  //DmObj *chain;

  pm_DbRecordIntGet (line, 9,  10, res_seq);
  pm_DbRecordIntGet (line, 14, 17, num_res);
  chain_id = line[11];

  //residues = dm_MemAlloc (0, PmResName, num_res);
  n = 0;

  while (1) {
    pm_DbRecordStrGet (line, 20, 70, str);
     s = (char *)strtok (str, " ");

    while (s && (*s != '\n') && (*s != ' ')) {
      //strcpy (residues[n++], s);
      s = (char *)strtok (NULL, " ");
      }

    if (n == num_res) {
      break;
      }

    if (!pm_DbLineGet(this->fp, line)) {
      break;
      }
    }
  }


//*============================================================*
//*==========           pm_DbPdbHelixProc            ==========*
//*============================================================*
// process HELIX entry

void
PmDbPdbInterface::pm_DbPdbHelixProc(char *line)
  {
  PmHelix helix;
  PmMolecule *mol;
  char id[4];

  /*
  fprintf (stderr, "\n ----------  pm_DbPdbHelixProc ---------- \n");
  fprintf (stderr, " line [%s] \n", line);
  */

  pm_DbRecordIntGet (line, 8, 10, helix.num);
  pm_DbRecordStrGet (line, 12, 14, id);
  pm_DbRecordStrGet (line, 16, 18, helix.init_res_name);
  pm_DbRecordStrGet (line, 20, 20, &helix.init_chain_id);
  pm_DbRecordIntGet (line, 22, 25, helix.init_seq_num);

  pm_DbRecordStrGet (line, 28, 30, helix.term_res_name);
  pm_DbRecordStrGet (line, 32, 32, &helix.term_chain_id);
  pm_DbRecordIntGet (line, 34, 37, helix.term_seq_num);

  pm_DbRecordIntGet (line, 39, 40, helix.hclass);
  pm_DbRecordIntGet (line, 72, 76, helix.length);

  for (unsigned int i = 0; i < strlen(id); i++) {
    if (id[i] != ' ') {
      helix.id.push_back(id[i]);  
      }
    }

  //fprintf (stderr, " >>> helix[%s] \n", helix.id.c_str());
  getCurrentMolecule(PM_MOLECULE_PROTEIN, &mol);
  mol->addHelix(helix);
  }

//*============================================================*
//*==========           pm_DbPdbHelixProc            ==========*
//*============================================================*
// process SHEET entry

void
PmDbPdbInterface::pm_DbPdbSheetProc(char *line)
  {
  PmSheet sheet;
  PmMolecule *mol;
  char id[4];

  /*
  fprintf (stderr, "\n ----------  pm_DbPdbSheetProc ---------- \n");
  fprintf (stderr, " line [%s] \n", line);
  */

  pm_DbRecordIntGet (line, 8, 10, sheet.num);
  pm_DbRecordStrGet (line, 12, 14, id);
  pm_DbRecordIntGet (line, 15, 16, sheet.num_strands);

  pm_DbRecordStrGet (line, 18, 20, sheet.init_res_name);
  pm_DbRecordStrGet (line, 22, 22, &sheet.init_chain_id);
  pm_DbRecordIntGet (line, 23, 26, sheet.init_seq_num);

  pm_DbRecordStrGet (line, 29, 31, sheet.term_res_name);
  pm_DbRecordStrGet (line, 33, 33, &sheet.term_chain_id);
  pm_DbRecordIntGet (line, 34, 37, sheet.term_seq_num);

  pm_DbRecordIntGet (line, 39, 40, sheet.sense);

  for (unsigned int i = 0; i < strlen(id); i++) {
    if (id[i] != ' ') {
      sheet.id.push_back(id[i]);
      }
    }

  //fprintf (stderr, " >>> add sheet[%s] \n", sheet.id.c_str());
  getCurrentMolecule(PM_MOLECULE_PROTEIN, &mol);
  mol->addSheet(sheet);
  }

