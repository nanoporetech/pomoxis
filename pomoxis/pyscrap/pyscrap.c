#include <Python.h>

#include <assert.h>
#include <dirent.h>
#include <err.h>
#include <glob.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>

#include "decode.h"
#include "networks.h"
#include "util.h"


static PyMethodDef module_functions[] = {
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(pyscrap)
{
    PyObject *m;

    MOD_DEF(m, "pyscrap", "Crappie scrappie wrappie",
            module_functions)

    if (m == NULL)
        return MOD_ERROR_VAL;
    return MOD_SUCCESS_VAL(m);

}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#  ifdef MODULE_API_EXPORTS
#    define MODULE_API __declspec(dllexport)
#    define restrict __restrict
#  else
#    define MODULE_API __declspec(dllimport)
#  endif
#else
#  define MODULE_API
#endif

MODULE_API int module_init();

#ifdef __cplusplus
}
#endif


MODULE_API void free_event_table(event_table et){
	free(et.event);
}

MODULE_API void free_char(char* ptr){
	free(ptr);
}

struct basecall_result {
	float score;
	char * bases;
};

MODULE_API struct basecall_result basecall_events(event_t * events, const int n_events){
        event_table et = (event_table){n_events, 0, n_events, events};

	const float min_prob = 1e-5;
	const float skip_pen = 0.0;
	const bool use_slip = false;
	const bool dwell_correction = false;

	scrappie_matrix post = nanonet_posterior(et, min_prob, true);
	if(NULL == post){
		printf("Post was NULL\n");
		return (struct basecall_result){0, NULL};
	}
	const int nev = post->nc;
	const int nstate = post->nr;

	int * history_state = calloc(nev, sizeof(int));
	float score = decode_transducer(post, skip_pen, history_state, use_slip);
	post = free_scrappie_matrix(post);
	int * pos = calloc(nev, sizeof(int));
	char * basecall = overlapper(history_state, nev, nstate - 1, pos);
	const size_t basecall_len = strlen(basecall);

	const int evoffset = et.start;
	for(int ev=0 ; ev < nev ; ev++){
		et.event[ev + evoffset].state = 1 + history_state[ev];
		et.event[ev + evoffset].pos = pos[ev];
	}

	if(dwell_correction){
		char * newbasecall = homopolymer_dwell_correction(et, history_state, nstate, basecall_len);
		if(NULL != newbasecall){
			free(basecall);
			basecall = newbasecall;
		}
	}
	free(pos);
	free(history_state);
	return (struct basecall_result){score, basecall};
}
