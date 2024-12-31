## Build

### MSM

To run the code:

```shell
cd demo/MSM
make thread
 ./build/test_pip_threads 
```

If there is a segment fault, then

```shell
ulimit -s "栈的大小"
```

