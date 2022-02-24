#ifndef C_STRING_HASH_EXTENSION
#define C_STRING_HASH_EXTENSION
struct cstrhash
{
    size_t operator()(const char *str) const
    {
        unsigned int hash = 0;
        unsigned int x = 0;
        while (*str)
        {
            hash = (hash << 4u) + (*str++);
            if ((x = hash & 0xF0000000u) != 0)
            {
                hash ^= (x >> 24u);
                hash &= ~x;
            }
        }
        return (hash & 0x7FFFFFFFu);
    }
};
struct cstrcmp {bool operator()(const char *s1, const char *s2) const {return strcmp(s1, s2) == 0;}};
#endif
