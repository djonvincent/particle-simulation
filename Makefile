CC=g++
CFLAGS=-O3 --std=c++11
TARGET=step4

SNAPSHOT=0.01
FINAL_TIME=0.05
DT=0.01
NUM_PARTICLES=10000

compile:
	$(CC) $(CFLAGS) solution-$(TARGET).c -o $(TARGET) $(OFLAGS)

.PHONY: run
run: compile genconds
	time ./$(TARGET) $(SNAPSHOT) $(FINAL_TIME) $(DT) initial-conditions.txt

genconds: 
	python create-initial-conditions.py $(NUM_PARTICLES) 1 1

.PHONY: clean
clean:
	rm $(TARGET)

