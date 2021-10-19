# Cube


Cube is a secondary development based on RTKLIB and mainly composed of two modules, satellite clock/decoupled clock estimation (the server end) and PPP/PPP-AR (the user end). At present, only GPS observations can be processed, and we only experimented with Cube on Windows system. However, benefiting from the C language, Cube can be easily transplanted to other operating systems, such as Linux and Macintosh. At the server end, Cube can estimate IGS legacy clocks, IGS legacy clocks with satellite code bias (SCB) extraction and decoupled clocks using undifferenced observations. Additionally, PPP-AR based on the IRC model and the DCK model can be implemented at the user end. Furthermore, the crucial parameters, such as receiver coordinates, the zenith tropospheric delays (ZTDs), ambiguities and residuals, are all output into formatted files, which is very convenient for additional applications and for post-analyses of the results.

## History
2021/10/15	1.0	new

## Contact

Shuai Liu

[lsnav@foxmail.com](mailto:lsnav@foxmail.com)


