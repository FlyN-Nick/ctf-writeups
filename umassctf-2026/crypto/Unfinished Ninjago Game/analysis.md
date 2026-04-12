### Overall Structure

The challenge first reads the flag, XORd with getrandom(&s,0x40,0).
It then gives you a menu option. You can to r to do a right thing, l to do a left. Each of these is basically the same, picks a random jump from PRNG.h and applying it to the output (after you send an A to apply) Then we get a leak of
s % 0x65 + DAT_0010d048 % 0x65 + DAT_0010d050 % 0x65 + DAT_0010d058 % 0x65 +
DAT_0010d060 % 0x65 + DAT_0010d068 % 0x65 + DAT_0010d070 % 0x65 + DAT_0010d078 % 0x65;
  putc((int)uVar1 + (int)(uVar1 / 0x65) * -0x65,stdout);

If we send an m

With these we should be able to fully constrain the initial state.

There is this other weird thing with the other inputs, you can read and wr

### Solution Idea
We can leak some state, then try to get the initial prng state out, which we can xor with the flag.

This should be doable since all the modifications are linear and thus simple to reverse in theory

We may also be able to do it with the poly arithmetic.
