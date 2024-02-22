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
      <anchor>aa5d2b27d9321575f1f6bce589aa0b2f9</anchor>
      <arglist>(int dtype, int64_t *stack, int64_t *heap)</arglist>
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
      <type>void</type>
      <name>PPIDD_Rank</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac713dfcfe82baefb567e15663332c98d</anchor>
      <arglist>(int64_t *me)</arglist>
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
      <anchor>ae15c21543c2b984588c32d32aff92f8b</anchor>
      <arglist>(void *buf, int64_t *count, int dtype, int64_t *dest, int64_t *sync)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Recv</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a30bb9805580edd9284e29313377b6bac</anchor>
      <arglist>(void *buf, int64_t *count, int dtype, int64_t *source, int64_t *lenreal, int64_t *sourcereal, int64_t *sync)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a6aa4dbb6fc408c0206af8a693dea1a20</anchor>
      <arglist>(int64_t *nodesel)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Iprobe</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8ab04d1b60d69a251463b12bfd7bde41</anchor>
      <arglist>(int64_t *tag, int64_t *source)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_BCast</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aef73e1faec9468f940eb3980753c9e76</anchor>
      <arglist>(void *buffer, int64_t *count, int dtype, int64_t *root)</arglist>
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
      <anchor>afea7848ae19f1a7942001f44fdc426f0</anchor>
      <arglist>(int dtype, void *buffer, int64_t *len, char *op)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Create_irreg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a1958bbfc0da771fbdca401e5edf69aab</anchor>
      <arglist>(char *name, int64_t *lenin, int64_t *nchunk, int dtype, int64_t *storetype, int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Create</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a42cbb824c3414ff606cb6f49e4bad003</anchor>
      <arglist>(char *name, int64_t *lentot, int dtype, int64_t *storetype, int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Destroy</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af8a751d551ce0e3913d2fe180170611e</anchor>
      <arglist>(int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Distrib</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a7b7cc054cc05957ef9a02ffa3910e6e8</anchor>
      <arglist>(int64_t *handle, int64_t *rank, int64_t *ilo, int64_t *ihi)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Location</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ae7f2178fab45fdcffb46c294010f4e6f</anchor>
      <arglist>(int64_t *handle, int64_t *ilo, int64_t *ihi, int64_t *map, int64_t *proclist, int64_t *np)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Get</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a3dc9a44d21d9cf560d19084919961eb1</anchor>
      <arglist>(int64_t *handle, int64_t *ilo, int64_t *ihi, void *buff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Put</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a0bcb9506432439ab6bd8b5aebb272f1e</anchor>
      <arglist>(int64_t *handle, int64_t *ilo, int64_t *ihi, void *buff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Acc</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aecaa019f7ed23c05a03c8b66a60fb473</anchor>
      <arglist>(int64_t *handle, int64_t *ilo, int64_t *ihi, void *buff, void *fac)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Read_inc</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a86cdc1a370b6ebd3b97fea4aac89f0cb</anchor>
      <arglist>(int64_t *handle, int64_t *inum, int64_t *incr, int64_t *returnval)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Zero_patch</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8f8345ae541d1fdcec2d1976cbf38267</anchor>
      <arglist>(int64_t *handle, int64_t *ilo, int64_t *ihi)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Zero</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aa71e4ea4a5318c48f037ed20097953b3</anchor>
      <arglist>(int64_t *handle)</arglist>
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
      <name>PPIDD_Duplicate</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a8fb85806da2adf2a61f57057fb4782f9</anchor>
      <arglist>(int64_t *handlei, int64_t *handlej, char *name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Inquire_name</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a797999e233bfa7facd26309045be55d5</anchor>
      <arglist>(int64_t *handle, char *name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Inquire_stype</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>abbf907f0689650965e18e1c731bd06c6</anchor>
      <arglist>(int64_t *handle, int64_t *storetype)</arglist>
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
      <anchor>a3fd148b9a02822e27567d2d5aa69e954</anchor>
      <arglist>(int64_t *storetype, int64_t *number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Lock_mutex</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a2b1aced8474d99c72ac36b5a818e9f7a</anchor>
      <arglist>(int64_t *inum)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Unlock_mutex</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a278aad6be9ad4e26bd741cf6ce43e0da</anchor>
      <arglist>(int64_t *inum)</arglist>
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
      <anchor>a1af52964bf4d63d1a1a99469784e9a23</anchor>
      <arglist>(char *name, int64_t *type, int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_write</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a1c63259ed65137b4e1a3c94ab220cf54</anchor>
      <arglist>(int64_t *handle, double *byte_offset, void *buff, int64_t *byte_length)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_awrite</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a00deb7f80cfb7672170a98340e4fb67a</anchor>
      <arglist>(int64_t *handle, double *byte_offset, void *buff, int64_t *byte_length, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_read</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac7fcd84c0746c5aedd42df88d361b8ca</anchor>
      <arglist>(int64_t *handle, double *byte_offset, void *buff, int64_t *byte_length)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_aread</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ad02d869076f1217b6d8aa68ced55db54</anchor>
      <arglist>(int64_t *handle, double *byte_offset, void *buff, int64_t *byte_length, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ad631bb122d34f165aee06a0efc4c805f</anchor>
      <arglist>(int64_t *handle, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_waitall</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a20ca98f4fc97804ed43deebc2cc9b1bd</anchor>
      <arglist>(int64_t *list, int64_t *num)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_probe</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aaccfa81c2ef7578941bcd838bfff84be</anchor>
      <arglist>(int64_t *request_id, int64_t *status)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_close</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ab325cb3fe584b7e16c265ece9575849e</anchor>
      <arglist>(int64_t *handle)</arglist>
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
      <anchor>aa4b693376a460bec0e0637b456157f8f</anchor>
      <arglist>(int64_t *handle, double *fsize)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Eaf_truncate</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aeac2bfe2f1c8c21ef63b7a52c92c7612</anchor>
      <arglist>(int64_t *handle, double *offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Eaf_errmsg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a58357b141947198c07d8133a46f19201</anchor>
      <arglist>(int *code, char *message)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_create</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>af35b2fa76ca10f0bfebaa096540a9b56</anchor>
      <arglist>(char *name, double *size_hard_limit, double *size_soft_limit, double *req_size, int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_write</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a3923188d3ab54e46d1eed4801b6f2135</anchor>
      <arglist>(int64_t *handle, double *byte_offset, double *byte_length, double *buff, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_read</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>a0ad5690a0214a5e22fa7211f17a46a43</anchor>
      <arglist>(int64_t *handle, double *byte_offset, double *byte_length, double *buff, int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_wait</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aca8c042520146598aaccc5b2a47a5dd8</anchor>
      <arglist>(int64_t *request_id)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_waitall</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>acb8ad3ce211ac24b6abe1fa7e959098d</anchor>
      <arglist>(int64_t *list, int64_t *num)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PPIDD_Sf_destroy</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>aa144a963d82f67586269a22ffd6432df</anchor>
      <arglist>(int64_t *handle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PPIDD_Sf_errmsg</name>
      <anchorfile>ppidd_8cpp.html</anchorfile>
      <anchor>ac288246d854eba6ded231d42528a6aa4</anchor>
      <arglist>(int *code, char *message)</arglist>
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
