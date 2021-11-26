
exec: main.c tasks.o
	gcc main.c tasks.o -o exec -lm -g
	gcc -c tasks.c

tasks.o: tasks.c
	gcc -c tasks.c

clean:
	rm *.o main