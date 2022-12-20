CC=gcc
BASEFLAGS= -Wall -lm -lpthread
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
CPRFFLAGS= $(BASEFLAGS) -pg
DEPS=truth_assign.h llist.h vector_ops.h lor.h ray_trace.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

containement: brightest_scatters.o llist.o vector_ops.o
	$(CC) -o containment/contained brightest_scatters.o llist.o vector_ops.o $(CFLAGS)

reverse_kin: inverse_kinematics.o llist.o vector_ops.o
	$(CC) -o reverse_kinematics $^ $(CFLAGS)

event_view: output_scatter_sets.o llist.o vector_ops.o
	$(CC) -o event_viewer $^ $(CFLAGS)

debug_reverse_kin: inverse_kinematics.c llist.c vector_ops.c
	$(CC) -o debug_reverse_kinematics $^ -Wall -lm -g

render: lor_render.o vector_ops.o ray_trace.o llist.o
	$(CC) -o renderer $^ $(CFLAGS)

traverse_render: lor_render_traversal.o vector_ops.o ray_trace.o llist.o
	$(CC) -o traverse_renderer $^ $(CFLAGS)

rod_traversal: rod_traversal.o vector_ops.o ray_trace.o llist.o
	$(CC) -o rod_traversal $^ $(CFLAGS)

light_reverse_kin: inverse_kinematics_lightweight.o llist.o vector_ops.o
	$(CC) -o light_reverse_kin $^ $(CFLAGS)

debug_light_reverse_kin: inverse_kinematics_lightweight.c llist.c vector_ops.c
	$(CC) -o debug_light_reverse_kinematics $^ $(CDBGFLAGS)

debug_render: lor_render.c vector_ops.c ray_trace.o llist.o
	$(CC) -o debug_renderer $^ $(CDBGFLAGS)

debug_traverse_render: lor_render_traversal.c vector_ops.c ray_trace.c llist.c
	$(CC) -o debug_traverse_renderer $^ $(CDBGFLAGS)

profile_light_reverse_kin: inverse_kinematics_lightweight.c llist.c vector_ops.c
	$(CC) -o profile_light_reverse_kinematics $^ $(CPRFFLAGS)

profile_reverse_kin: inverse_kinematics.c llist.c vector_ops.c
	$(CC) -o profile_reverse_kinematics $^ $(CPRFFLAGS)

debug_ray_trace: ray_trace.c ray_test.c vector_ops.c
	$(CC) -o debug_ray_trace $^ $(CDBGFLAGS)

classical_pet: classical_pet.o llist.o vector_ops.o
	$(CC) -o classical_pet $^ $(CFLAGS)

debug_classical_pet: classical_pet.c llist.c vector_ops.c
	$(CC) -o debug_classical_pet $^ $(CDBGFLAGS)