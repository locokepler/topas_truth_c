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
	$(CC) -o reverse_kinimatics $^ $(CFLAGS)

event_view: output_scatter_sets.o llist.o vector_ops.o
	$(CC) -o event_viewer $^ $(CFLAGS)

debug_reverse_kin: inverse_kinematics.c llist.c vector_ops.c
	$(CC) -o debug_reverse_kinimatics $^ -Wall -lm -g

render: lor_render.o vector_ops.o
	$(CC) -o renderer $^ $(CFLAGS)

light_reverse_kin: inverse_kinematics_lightweight.o llist.o vector_ops.o
	$(CC) -o light_reverse_kin $^ $(CFLAGS)

debug_light_reverse_kin: inverse_kinematics_lightweight.c llist.c vector_ops.c
	$(CC) -o debug_light_reverse_kinematics $^ $(CDBGFLAGS)

debug_render: lor_render.c vector_ops.c
	$(CC) -o debug_renderer $^ $(CDBGFLAGS)

profile_light_reverse_kin: inverse_kinematics_lightweight.c llist.c vector_ops.c
	$(CC) -o profile_light_reverse_kinematics $^ $(CPRFFLAGS)

debug_ray_trace: ray_trace.c ray_test.c vector_ops.c
	$(CC) -o debug_ray_trace $^ $(CDBGFLAGS)