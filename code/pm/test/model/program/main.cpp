
#include "pm.h"
#include "db.h"
#include "mol.h"
#include "rbsim.h"

using namespace ProteinModeler;


int
main (int nargs, char **args) 
  {


  int num_coords;

  PmDbInterfaceSelect db_sel;

  PmDbInterface *db;

  PmAtoms *atoms;

  PmMolecule *mol;

  PmRigidSimulation *sim;


  // read molecule // 

  PmDbType type;

  type = pm_DbFormatConv ("pdb");
  char *name = "iq.pdb";
  db = db_sel.create (type); 
  db->open (name);
  db->read ("mol1");

  fprintf (stderr, "\n"); 
  mol = db->getMolecule();
  fprintf (stderr, ">>>>>> num mol atoms = [%d] \n", mol->getNumAtoms());

  int num_res = mol->getNumResidues ();
  fprintf (stderr, ">>>>>> num mol res [%d] \n", num_res); 


  // create domains //

  name = "d1";
  char *desc1 = "A";
  PmMolecule *domainA = mol->createDomain (name, desc1); 

  name = "d2";
  char *desc2 = "B";
  PmMolecule *domainB = mol->createDomain (name, desc2); 


  // create a rigid body simulation //

  sim = new PmRigidSimulation ("rigid sim");

  // create and add two rigid bodies
  PmRigidBody *body1 = new PmRigidBody ("body1", domainA); 
  sim->addBody (body1);

  PmRigidBody *body2 = new PmRigidBody ("body2", domainB); 
  sim->addBody (body2);

  // create and add a joint connecting the two rigid bodies
  PmJoint *joint1 = new PmJoint ("body1-body2-hinge", PM_JOINT_PIN); 
  PmVector3 pos;
  pos.set (1,1,1);
  joint1->setPosition (pos);
  joint1->setBodies (body1, body2);
  sim->addJoint (joint1);


  }


