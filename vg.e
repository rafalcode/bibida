==27528== Memcheck, a memory error detector
==27528== Copyright (C) 2002-2013, and GNU GPL'd, by Julian Seward et al.
==27528== Using Valgrind-3.10.0.SVN and LibVEX; rerun with -h for copyright info
==27528== Command: /home/nutria/mylocal/bin/cdsck2_d cele_x.fa cel.fa
==27528== 
--27528-- Valgrind options:
--27528--    --leak-check=full
--27528--    --show-reachable=yes
--27528--    --track-origins=yes
--27528--    -v
--27528-- Contents of /proc/version:
--27528--   Linux version 3.13.0-49-generic (buildd@akateko) (gcc version 4.8.2 (Ubuntu 4.8.2-19ubuntu1) ) #83-Ubuntu SMP Fri Apr 10 20:11:33 UTC 2015
--27528-- Arch and hwcaps: AMD64, amd64-cx16-lzcnt-rdtscp-sse3-avx-avx2-bmi
--27528-- Page sizes: currently 4096, max supported 4096
--27528-- Valgrind library directory: /usr/lib/valgrind
--27528-- Reading syms from /home/nutria/mylocal/bin/cdsck2_d
--27528-- Reading syms from /lib/x86_64-linux-gnu/ld-2.19.so
--27528--   Considering /lib/x86_64-linux-gnu/ld-2.19.so ..
--27528--   .. CRC mismatch (computed 4cbae35e wanted 8d683c31)
--27528--   Considering /usr/lib/debug/lib/x86_64-linux-gnu/ld-2.19.so ..
--27528--   .. CRC is valid
--27528-- Reading syms from /usr/lib/valgrind/memcheck-amd64-linux
--27528--   Considering /usr/lib/valgrind/memcheck-amd64-linux ..
--27528--   .. CRC mismatch (computed 37cdde19 wanted adc367dd)
--27528--    object doesn't have a symbol table
--27528--    object doesn't have a dynamic symbol table
--27528-- Scheduler: using generic scheduler lock implementation.
--27528-- Reading suppressions file: /usr/lib/valgrind/default.supp
==27528== embedded gdbserver: reading from /tmp/vgdb-pipe-from-vgdb-to-27528-by-nutria-on-???
==27528== embedded gdbserver: writing to   /tmp/vgdb-pipe-to-vgdb-from-27528-by-nutria-on-???
==27528== embedded gdbserver: shared mem   /tmp/vgdb-pipe-shared-mem-vgdb-27528-by-nutria-on-???
==27528== 
==27528== TO CONTROL THIS PROCESS USING vgdb (which you probably
==27528== don't want to do, unless you know exactly what you're doing,
==27528== or are doing some strange experiment):
==27528==   /usr/lib/valgrind/../../bin/vgdb --pid=27528 ...command...
==27528== 
==27528== TO DEBUG THIS PROCESS USING GDB: start GDB like this
==27528==   /path/to/gdb /home/nutria/mylocal/bin/cdsck2_d
==27528== and then give GDB the following command
==27528==   target remote | /usr/lib/valgrind/../../bin/vgdb --pid=27528
==27528== --pid is optional if only one valgrind process is running
==27528== 
--27528-- REDIR: 0x4019ca0 (strlen) redirected to 0x38068331 (???)
--27528-- Reading syms from /usr/lib/valgrind/vgpreload_core-amd64-linux.so
--27528--   Considering /usr/lib/valgrind/vgpreload_core-amd64-linux.so ..
--27528--   .. CRC mismatch (computed 329d6860 wanted c0186920)
--27528--    object doesn't have a symbol table
--27528-- Reading syms from /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so
--27528--   Considering /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so ..
--27528--   .. CRC mismatch (computed 1fb85af8 wanted 2e9e3c16)
--27528--    object doesn't have a symbol table
==27528== WARNING: new redirection conflicts with existing -- ignoring it
--27528--     old: 0x04019ca0 (strlen              ) R-> (0000.0) 0x38068331 ???
--27528--     new: 0x04019ca0 (strlen              ) R-> (2007.0) 0x04c2e1a0 strlen
--27528-- REDIR: 0x4019a50 (index) redirected to 0x4c2dd50 (index)
--27528-- REDIR: 0x4019c70 (strcmp) redirected to 0x4c2f2f0 (strcmp)
--27528-- REDIR: 0x401a9c0 (mempcpy) redirected to 0x4c31da0 (mempcpy)
--27528-- Reading syms from /lib/x86_64-linux-gnu/libc-2.19.so
--27528--   Considering /lib/x86_64-linux-gnu/libc-2.19.so ..
--27528--   .. CRC mismatch (computed dc620abc wanted 148cbd6e)
--27528--   Considering /usr/lib/debug/lib/x86_64-linux-gnu/libc-2.19.so ..
--27528--   .. CRC is valid
--27528-- REDIR: 0x4ec3d60 (strcasecmp) redirected to 0x4a25720 (_vgnU_ifunc_wrapper)
--27528-- REDIR: 0x4ec6050 (strncasecmp) redirected to 0x4a25720 (_vgnU_ifunc_wrapper)
--27528-- REDIR: 0x4ec3530 (memcpy@GLIBC_2.2.5) redirected to 0x4a25720 (_vgnU_ifunc_wrapper)
--27528-- REDIR: 0x4ec17c0 (rindex) redirected to 0x4c2da30 (rindex)
--27528-- REDIR: 0x4ebe070 (strcmp) redirected to 0x4a25720 (_vgnU_ifunc_wrapper)
--27528-- REDIR: 0x4f77200 (__strcmp_ssse3) redirected to 0x4c2f1b0 (strcmp)
--27528-- REDIR: 0x4eba220 (calloc) redirected to 0x4c2cbf0 (calloc)
--27528-- REDIR: 0x4eb9750 (malloc) redirected to 0x4c2ab10 (malloc)
--27528-- REDIR: 0x4ec2410 (__GI_strstr) redirected to 0x4c32030 (__strstr_sse2)
--27528-- REDIR: 0x4eb9ef0 (realloc) redirected to 0x4c2ce10 (realloc)
--27528-- REDIR: 0x4ec35c0 (memset) redirected to 0x4c31350 (memset)
--27528-- REDIR: 0x4ebfee0 (strncmp) redirected to 0x4a25720 (_vgnU_ifunc_wrapper)
--27528-- REDIR: 0x4f78460 (__strncmp_ssse3) redirected to 0x4c2e8c0 (strncmp)
--27528-- REDIR: 0x4eb9df0 (free) redirected to 0x4c2bd80 (free)
--27528-- REDIR: 0x4ecaac0 (strchrnul) redirected to 0x4c319b0 (strchrnul)
--27528-- REDIR: 0x4ec8780 (__GI_memcpy) redirected to 0x4c2fc90 (__GI_memcpy)
--27528-- REDIR: 0x4ebfac0 (strlen) redirected to 0x4c2e0e0 (strlen)
==27528== Invalid free() / delete / delete[] / realloc()
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403EBA: main (cdsck2.c:571)
==27528==  Address 0x5298d70 is 0 bytes inside a block of size 24 free'd
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403C08: main (cdsck2.c:553)
==27528== 
==27528== Invalid free() / delete / delete[] / realloc()
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403EED: main (cdsck2.c:572)
==27528==  Address 0x52df910 is 0 bytes inside a block of size 1,456 free'd
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403C3B: main (cdsck2.c:554)
==27528== 
cele_x.fa/numseq=30/numfiltered=26/%filt=86.7
cel.fa/numseq=3/numfiltered=0/%filt=0.0
==27528== 
==27528== HEAP SUMMARY:
==27528==     in use at exit: 0 bytes in 0 blocks
==27528==   total heap usage: 14,236 allocs, 14,296 frees, 24,466,148 bytes allocated
==27528== 
==27528== All heap blocks were freed -- no leaks are possible
==27528== 
==27528== ERROR SUMMARY: 60 errors from 2 contexts (suppressed: 0 from 0)
==27528== 
==27528== 30 errors in context 1 of 2:
==27528== Invalid free() / delete / delete[] / realloc()
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403EED: main (cdsck2.c:572)
==27528==  Address 0x52df910 is 0 bytes inside a block of size 1,456 free'd
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403C3B: main (cdsck2.c:554)
==27528== 
==27528== 
==27528== 30 errors in context 2 of 2:
==27528== Invalid free() / delete / delete[] / realloc()
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403EBA: main (cdsck2.c:571)
==27528==  Address 0x5298d70 is 0 bytes inside a block of size 24 free'd
==27528==    at 0x4C2BDEC: free (in /usr/lib/valgrind/vgpreload_memcheck-amd64-linux.so)
==27528==    by 0x403C08: main (cdsck2.c:553)
==27528== 
==27528== ERROR SUMMARY: 60 errors from 2 contexts (suppressed: 0 from 0)
