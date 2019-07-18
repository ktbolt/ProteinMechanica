
src = \
      atom.cpp atoms.cpp \
      bexp.cpp bez.cpp body.cpp bsim.cpp  \
      cmpd.cpp cmd.cpp csys.cpp \
      db.cpp  \
      force.cpp \
      gc.cpp geom.cpp graphics.cpp grid.cpp \
      joint.cpp     \
      mesh.cpp model.cpp mol.cpp motor.cpp msg.cpp msr.cpp  \
      ode_solv.cpp              \
      part.cpp pm.cpp pobj.cpp pot.cpp \
      rbsim.cpp rbsolv.cpp res.cpp rest.cpp    \
      sd.cpp sim.cpp sobj.cpp solid.cpp state.cpp surf.cpp pmsys.cpp \
      trace.cpp

obj = $(src:.cpp=.o)
cmd_deps = cmd_body.cpp cmd_bsim.cpp cmd_db.cpp cmd_dom.cpp cmd_force.cpp cmd_gr.cpp \
           cmd_curve.cpp cmd_grid.cpp    \
           cmd_joint.cpp cmd_mbody.cpp cmd_model.cpp cmd_motor.cpp cmd_mol.cpp          \
           cmd_msr.cpp cmd_part.cpp cmd_pot.cpp cmd_sim.cpp cmd_solid.cpp cmd_surf.cpp \
           cmd_sys.cpp cmd_trace.cpp cmd_units.cpp cmd_vector.cpp



