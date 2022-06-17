CC=gcc
CFLAGS=-Wall -lm -O2 -lpthread
DEPS=truth_assign.h llist.h vector_ops.h lor.h

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

debug_render: lor_render.c vector_ops.c
	$(CC) -o debug_renderer $^ -Wall -lm -lpthread -g