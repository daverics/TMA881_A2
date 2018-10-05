char *concatDiffString(int argc, char *argv[]) {
     int length = 0;
     for (int i = 0; i < argc; ++i)
         length += strlen(argv[i]);


     char *output = (char*) malloc(length + 1);


     char *dest = output;
     for (int i = 0; i < argc; ++i) {
         char *src = argv[i];
         while (*src)
             *dest++ = *src++;
     }
     *dest = '\0';
     return output;
}
