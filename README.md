# Go-Smart Goosefoot Control Tool

This repository contains a series of helper scripts to support users examining output of containerized Goosefoot simulations.
This is intended to be cloned behind-the-scenes by [glot](https://github.com/go-smart/glot), a tool that manages a local Glossia
server, as well as aiding in remote analysis.

## Re-running

One of the key features of the Glossia system is the ability to run locally, retrospectively your simulation in an identical
environment and configuration to that executed on a remote Glossia server. This works by recording the container image ID used
(`input\glossia-simimage.txt`) in the input directory. When a diagnostic bundle is requested from Glossia, all input to the
container is included, and the scripts in this repository give the user the option of running with the exact same simulation
image or the latest container image on their machine.

These instructions create necessary directories and execute the Goosefoot solver in the same container that is used on the Go-Smart server.

```bash
sudo ./setup.sh
sudo ./run.sh
```

Output appears in the `output` folder, created by `setup.sh`. Note that, as the container will create files as its own unprivileged user,
you may need to `sudo chown -R $(whoami) .` within the current directory.


## Workflow

Normal usage, including the wider tool-chain. `$SIMPREFIX` is a least the start of the simulation UUID of the simulation on the server.
The UUID may be discovered from, for instance, the output (on the server) of `glot search` or Glossia logs. If a third-party front-end is being used, it may support
user-friendly UUID reporting. The full UUID, which is `$SIMPREFIX` or, if `$SIMPREFIX` is less than 32 characters, has `$SIMPREFIX` as a prefix,
is denoted `$SIMUUID`.

### Checking results

* If Docker host is local (`inspect-diagnostic` optional)
** `glot results --inspect-diagnostic $SIMPREFIX`
** `cd $SIMUUID/output`
* If Docker host is remote
** `(remote) glot results $SIMPREFIX`
** Transmit `$SIMPREFIX-results.tgz` to user

### Debugging

* If Docker host is local
** `glot diagnostic --inspect $SIMPREFIX`
* If Docker host is remote
** `(remote) glot diagnostic $SIMPREFIX`
** Transmit `$SIMPREFIX-diagnostic.tgz` to user
** `(local) glot inspect $SIMPREFIX-diagnostic.tgz`
* `cd $SIMUUID`
* Check `output/logs` for information
* `sudo ./setup.sh && sudo ./run.sh`
* Confirm output is as expected
* To do non-containerized debugging afterwards, `cd output/run && chown -R $(whoami) .` will give you a normal (used) simulation working directory.

Note that results provided from the Glossia are only those sent back from the container, so re-running will be necessary to obtain intermediary
data.
