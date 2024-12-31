CFLAGS=-O3 -march=native -mtune=native -g
LIBS=target/lib/librelic_s.a -lgmp
RELIC_ROOT = ../..

INC_DIR = include
SRC_DIR = src
BUILD_DIR = build
TEST_DIR = test

CDEFINES ?=

pip: 
	gcc ${CFLAGS} $(SRC_DIR)/data_struct.c $(SRC_DIR)/pippenger.c  $(TEST_DIR)/test_pippenger.c -o $(BUILD_DIR)/test_pippenger -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} ${LIBS} -lm

ifma:  
	gcc ${CFLAGS} $(SRC_DIR)/data_struct.c $(SRC_DIR)/gfparith.c $(SRC_DIR)/gfp2arith.c $(SRC_DIR)/padd.c $(SRC_DIR)/utils.c $(SRC_DIR)/pippenger.c $(SRC_DIR)/pip_ifma.c $(TEST_DIR)/test_pip_ifma.c -o $(BUILD_DIR)/test_pip_ifma -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} ${LIBS} -lm -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -fwrapv	$(CDEFINES)

pair_ifma: 
	gcc ${CFLAGS} $(SRC_DIR)/* $(TEST_DIR)/test_pair_ifma.c -o $(BUILD_DIR)/test_pair_ifma -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} ${LIBS} -lm -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -fwrapv $(CDEFINES)

thread: 
	gcc ${CFLAGS} $(SRC_DIR)/* $(TEST_DIR)/test_pip_threads.c -o $(BUILD_DIR)/test_pip_threads -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} ${LIBS} -lm -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -fwrapv -fopenmp $(CDEFINES)

pair_thread: 
	gcc ${CFLAGS} $(SRC_DIR)/* $(TEST_DIR)/test_pair_threads.c -o $(BUILD_DIR)/test_pair_threads -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} ${LIBS} -lm -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -fwrapv -fopenmp $(CDEFINES)

msm: libmsm.a

libmsm.a: $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
	ar rcs $(BUILD_DIR)/libmsm.a $(BUILD_DIR)/*.o target/lib/librelic_s.a
	rm -rf $(BUILD_DIR)/*.o

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	gcc -c ${CFLAGS} -I ${RELIC_ROOT}/include -I target/include -I ${INC_DIR} $< -o $@ -lgmp -mavx512ifma -fomit-frame-pointer -m64 -mbmi2 -fwrapv -fopenmp


lib:
	mkdir -p target
	cd target; ${RELIC_ROOT}/../preset/x64-pbc-bls12-381.sh ${RELIC_ROOT}/../; cmake -DEP_METHD='JACOB;LWNAF;COMBS;INTER;SSWUM' .; make

clean:
	rm -rf $(BUILD_DIR)/*


