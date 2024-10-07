<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.1">
  <compound kind="file">
    <name>ppidd.cpp</name>
    <path>/__w/ppidd/ppidd/src/</path>
    <filename>ppidd_8cpp.html</filename>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Initialize</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ad0bac13ed66299b75b584da2c0963ee0</anchor>
      <arglist>(int *argc, char ***argv, int impl)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Initialize_data</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>accd1aaec40977ac816e292441b1fdb6f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int64_t</type>
      <name>PPIDD_Worker_comm</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ab8a9bad19a3488c3c2eb193068dee536</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Finalize</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a768ee79b70214209cec6ad98fb0acdd6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Uses_ma</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a2a9621816492592025267948fc348b71</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_MA_init</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8a29bf48a081ca9d09da5cd415eda21f</anchor>
      <arglist>(int dtype, int64_t stack, int64_t heap)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PPIDD_Wtime</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8737999edc5cb2e59357e005a7d52c88</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Error</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a16e1aebb607213730b037da415349eba</anchor>
      <arglist>(char *message, int code)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Helper_server</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>afb8ee257625e38c67161287c34d5a028</anchor>
      <arglist>(int flag, int numprocs_per_server)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Size_all</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac1f35e2a2c40396f0de4be4c9894327f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Size</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>abe8f2ab3c7967ce078813ad0821e156f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Rank</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a26d9775043f3893d1ae8cb33f56912c7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Init_fence</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a5c67417ceed19a2f9e6bd5eb6bada439</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Fence</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a2aaef4140504ad475c1621e80ec3c774</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Send</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a32bf1e3984512dab4f53532232779675</anchor>
      <arglist>(void *buf, int count, int dtype, int dest, int sync)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Recv</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ae070a87257dfcf92320bea61fe50f0dd</anchor>
      <arglist>(void *buf, int count, int dtype, int source, int64_t *lenreal, int64_t *sourcereal, int sync)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a7d72d5ccab8bc41ee2e530a56a54fee3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Iprobe</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>abb13c8254fc47245acd9cdc8bedf08d2</anchor>
      <arglist>(int tag, int source)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_BCast</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a33e604b3772bc88bb0c41b8cd971142c</anchor>
      <arglist>(void *buffer, int count, int dtype, int root)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Barrier</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>abe47ed108d82f5d212dbd60cd13966b6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Gsum</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aa72378b9de01b4cae8ee2713549f1833</anchor>
      <arglist>(int dtype, void *buffer, int len, char *op)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Create_irreg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aabaa345ef135d2d65f32cb2de7461faa</anchor>
      <arglist>(char *name, int64_t *lenin, int64_t nchunk, int dtype, int storetype, int *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Create</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a0cf53817e9e57b8b3f807c2ca9ff613e</anchor>
      <arglist>(char *name, int64_t lentot, int dtype, int storetype, int *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Destroy</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af6eb8590dea63b899c9bca929848de80</anchor>
      <arglist>(int handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Distrib</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ad4f931760b39afbc5412f793bd337b28</anchor>
      <arglist>(int handle, int rank, int64_t *ilo, int64_t *ihi)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Location</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>afdcaa0a58adfdf57f71887d94704a756</anchor>
      <arglist>(int handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int *np)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Get</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a130e24887f4d0b506e7a5eda48d94ecd</anchor>
      <arglist>(int handle, int64_t *ilo, int64_t *ihi, void *buff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Put</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a9577128a5615c4776c76641a03293b74</anchor>
      <arglist>(int handle, int64_t *ilo, int64_t *ihi, void *buff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Acc</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac9efbd27a6608ae9df73f643f0400654</anchor>
      <arglist>(int handle, int64_t *ilo, int64_t *ihi, void *buff, void *fac)</arglist>
    </member>
    <member kind="function">
      <type>int64_t</type>
      <name>PPIDD_Read_inc</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a704ef6d0582d73b819007418973dea96</anchor>
      <arglist>(int handle, int64_t inum, int64_t incr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Zero_patch</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ab846031b89644e69555da8965f355f3f</anchor>
      <arglist>(int handle, int64_t ilo, int64_t ihi)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Zero</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a394fbda72a12f97a2cbeff94b79f484e</anchor>
      <arglist>(int handle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Nxtval</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a7f68db9670b16b6aa6bc2450f980a8ac</anchor>
      <arglist>(int numproc, int64_t *val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Inquire_name</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a704b674ff7f5dcc6b88c666c3308adef</anchor>
      <arglist>(int handle, char *name)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Inquire_stype</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a4e90163f6ad9674cd2086bd6c76e550e</anchor>
      <arglist>(int handle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Inquire_mem</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a6a56a2ff876976509db45c72aef808b1</anchor>
      <arglist>(int64_t *mem_used)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Create_mutexes</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a6c14480c1fe6301e5b53568cc82ccf75</anchor>
      <arglist>(int storetype, int number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Lock_mutex</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a36fddef2d864ffe78258c7ed40c2bf0f</anchor>
      <arglist>(int inum)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Unlock_mutex</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8c1b2bcce35bb241c10f6d9f470bdb7c</anchor>
      <arglist>(int inum)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Destroy_mutexes</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af4ee51d8abe235091c9a59863d466c80</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_open</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a88aa59519cbf4cb7b837d37be740af85</anchor>
      <arglist>(char *name, int type, int *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_write</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac778d6d42d0a1a06f92385eecbe6b8c0</anchor>
      <arglist>(int handle, double *byte_offset, void *buff, int64_t *byte_length)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_awrite</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a5aeae477abd2525befc0c158cb1e8702</anchor>
      <arglist>(int handle, double *byte_offset, void *buff, int64_t *byte_length, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_read</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a289185e09b3ac72eeb7c30839d05f57d</anchor>
      <arglist>(int handle, double *byte_offset, void *buff, int64_t *byte_length)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_aread</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a845d38877aea6a944107978f63fe53cf</anchor>
      <arglist>(int handle, double *byte_offset, void *buff, int64_t *byte_length, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a00e96f4c49ccd6ff69fa996cd26f8147</anchor>
      <arglist>(int handle, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_waitall</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aa006cca817f50c51c195285642e89ff9</anchor>
      <arglist>(int64_t *list, int num)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_probe</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af5f94372c339f1e4077ba41b517cf42d</anchor>
      <arglist>(int64_t *request_id, int *status)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_close</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a49b5d086751d2081fb99383b5ca3e130</anchor>
      <arglist>(int handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_delete</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a307d1d9646109797ab3fa84858012df4</anchor>
      <arglist>(char *name)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_length</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>add796b57fa6d7dc8862b22a900585a0a</anchor>
      <arglist>(int handle, double *fsize)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_truncate</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af132f8e9bc8bcd4aebad805a38477d35</anchor>
      <arglist>(int handle, double *offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Eaf_errmsg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>adf61b9e586d04071e76898a34f4cfca5</anchor>
      <arglist>(int code, char *message)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_create</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aa20fccc5bd21559a1423ca90d099f1f2</anchor>
      <arglist>(char *name, double size_hard_limit, double size_soft_limit, double req_size, int *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_write</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a01f04a1cb007eefeebc7a709b5e2c9dd</anchor>
      <arglist>(int handle, double byte_offset, double byte_length, double *buff, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_read</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a504ffa4bc617ba8c3c58e037fd7fbc82</anchor>
      <arglist>(int handle, double byte_offset, double byte_length, double *buff, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a94a234cd54f9d3575408f74fd6671d17</anchor>
      <arglist>(int64_t request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_waitall</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a17eaf202e52cc22757c6d11725b506ac</anchor>
      <arglist>(int64_t *list, int num)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_destroy</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ae2d5c28fddfa367c29ef79b4c9f3426f</anchor>
      <arglist>(int handle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Sf_errmsg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a9dd05422e0e5ad25f0b535e30616742b</anchor>
      <arglist>(int code, char *message)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ppidd_doxygen.h</name>
    <path>/__w/ppidd/ppidd/src/</path>
    <filename>ppidd__doxygen_8h.html</filename>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>PPIDD Reference Manual</title>
    <filename>index.html</filename>
  </compound>
</tagfile>
