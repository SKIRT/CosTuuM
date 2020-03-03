add_library(${LIBQUICKSCHED} STATIC qsched.c queue.c)
add_library(${LIBQUICKSCHED_SHARED} SHARED qsched.c queue.c)
