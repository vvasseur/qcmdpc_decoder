#define LOG2(v)                                                                \
    ((v) == 0                                                                  \
         ? 0                                                                   \
         : ((v) < 2                                                            \
                ? 1                                                            \
                : ((v) < 4                                                     \
                       ? 2                                                     \
                       : ((v) < 8                                              \
                              ? 3                                              \
                              : ((v) < 16                                      \
                                     ? 4                                       \
                                     : ((v) < 32                               \
                                            ? 5                                \
                                            : ((v) < 64                        \
                                                   ? 6                         \
                                                   : ((v) < 128                \
                                                          ? 7                  \
                                                          : ((v) < 256         \
                                                                 ? 8           \
                                                                 : 9)))))))))