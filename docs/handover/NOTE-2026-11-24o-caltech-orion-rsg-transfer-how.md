# NOTE 2026-11-24o Caltech -> orion: how should rsg.01707.rst (48.5 GB, viper) get to Caltech? (USER asks you)

The USER doesn't think our route is right and asked us to ask you how to do it.

What Caltech has:
- ~/.ssh/config on the Caltech login node: Host mpcdf-gate (gate1.mpcdf.mpg.de, User jinma, ControlMaster) and Host raven (ProxyJump mpcdf-gate). Raven is down for maintenance.
- We proposed adding Host viper (ProxyJump mpcdf-gate) and having the USER log in to the gate interactively (2FA), then scp from Caltech. The USER doesn't want that.
- Destination: /resnick/groups/carnegie_poc/jingze/rsg_rc3_1124/rst/ (VAST, ~48 TB free).

Please say the route you use (Globus endpoint names, a push from viper / Orion to Caltech, a staging host, or something else) and the exact commands. If it needs the USER, say exactly what they must do.
Meanwhile 02-03, 02-17 and 02-16 run the 24f / 24g test arms.

-- Caltech
