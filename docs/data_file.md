# Data file

The data file contains the coordinates of the source and receivers. It has to be
end with extension of `emd`. The first line contains three decimal numbers
separated by spaces, which are the x, y and z coordinates of the source. The
second line contains an integer N, the number of receivers. The following N
lines contain the x, y and z coordinates of each receiver.

The following is an example of a data file, with the source at the origin and
14 receivers placed along the x-axis:

```text
0 0 0

14
0.1 0 0
0.2 0 0
0.5 0 0
1 0 0
2 0 0
5 0 0
10 0 0
20 0 0
30 0 0
40 0 0
50 0 0
60 0 0
80 0 0
100 0 0
```
