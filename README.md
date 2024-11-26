To run:

```sh
docker build . -t spl-artifact
docker run -it spl-artifact:latest
```

Inside the docker image:
```sh
# setup environment
./start.sh
# run the build and benchmark
./run.sh
```
