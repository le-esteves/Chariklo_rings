export OPENMP=1
export PROFILING=0
export OPENGL=0# Set this to 1 to enable OpenGL
export SERVERL=0# Set this to 1 to enable the visualization web server
include ../../rebound_gods/src/Makefile.defs

RESTART_FLAG = 
LOG_FLAG =

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
SRC ?= problem.c

all: $(SRC) librebound
	$(CC) -I../../rebound_gods/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(SRC) -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful. To run REBOUND GODS VERSION, execute the file '$(EXEREBOUND)'."
	@mkdir -p Outputs

run:
	@./rebound
	
restart: RESTART_FLAG = -DRESTART
restart: $(SRC) librebound
	$(CC) -I../../rebound_gods/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) $(SRC) -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful using RESTARTFLAG. To run REBOUND GODS VERSION, execute the file '$(EXEREBOUND)'."

log: LOG_FLAG = -DLOG
log: $(SRC) librebound
	$(CC) -I../../rebound_gods/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) $(LOG_FLAG) $(SRC) -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful with LOGGING enabled. To run REBOUND GODS VERSION, execute the file '$(EXEREBOUND)'."
	@mkdir -p Outputs

restart_log: RESTART_FLAG = -DRESTART
restart_log: LOG_FLAG = -DLOG
restart_log: $(SRC) librebound
	$(CC) -I../../rebound_gods/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) $(LOG_FLAG) $(SRC) -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful with RESTART and LOGGING enabled. To run REBOUND GODS VERSION, execute the file '$(EXEREBOUND)'."

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../rebound_gods/src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@-$(RM) -r Outputs
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../rebound_gods/src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)

rebound_webgl.html: $(SRC)
	@echo "Compiling $(SRC) with emscripten (WebGL enabled)..."
	emcc -O3 -I../../rebound_gods/src/ ../../rebound_gods/src/*.c $(SRC) -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../rebound_gods/web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: $(SRC)
	@echo "Compiling $(SRC) with emscripten (WebGL disabled)..."
	emcc -O3 -I../../rebound_gods/src/ ../../rebound_gods/src/*.c $(SRC) -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../rebound_gods/web_client/shell_rebound_console.html -o rebound_console.html
