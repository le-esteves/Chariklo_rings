export OPENMP=1
export OPENGL=0# Set this to 1 to enable OpenGL
export SERVERL=1# Set this to 1 to enable the visualization web server
include ../../rebound_tp/src/Makefile.defs

RESTART_FLAG = 
LOG_FLAG =

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
all: problem.c librebound
	$(CC) -I../../rebound_tp/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@mkdir Outputs

run:
	@./rebound
	
restart: RESTART_FLAG = -DRESTART
restart: problem.c librebound
	$(CC) -I../../rebound_tp/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) problem.c -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful using RESTARTFLAG. To run REBOUND, execute the file '$(EXEREBOUND)'."

log: LOG_FLAG = -DLOG
log: problem.c librebound
	$(CC) -I../../rebound_tp/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) $(LOG_FLAG) problem.c -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful with LOGGING enabled. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@mkdir Outputs

restart_log: RESTART_FLAG = -DRESTART
restart_log: LOG_FLAG = -DLOG
restart_log: problem.c librebound
	$(CC) -I../../rebound_tp/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(RESTART_FLAG) $(LOG_FLAG) problem.c -L. -lrebound $(LIB) -o rebound
	@echo "Compilation successful with RESTART and LOGGING enabled. To run REBOUND, execute the file '$(EXEREBOUND)'."

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../rebound_tp/src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@-$(RM) -r Outputs
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../rebound_tp/src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../rebound_tp/src/ ../../rebound_tp/src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../rebound_tp/web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../rebound_tp/src/ ../../rebound_tp/src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../rebound_tp/web_client/shell_rebound_console.html -o rebound_console.html
