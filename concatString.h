char* concatString(int numDone, int numWrite) {
    int length = numWrite*strlen(colours[roots[0]]);

    char *output = (char*)malloc(length + 1);

    char *dest = output;
    for (int i = numDone; i < numDone + numWrite; ++i) {
         char *src = colours[roots[i]];
               while (*src)
               *dest++ = *src++;
               }
     *dest = '\0';
     return output;
}
