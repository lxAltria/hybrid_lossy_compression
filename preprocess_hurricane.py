import numpy as np

variables = np.array(["QCLOUDf48_log10", "QGRAUPf48_log10", "QICEf48_log10", "QRAINf48_log10", "QSNOWf48_log10", "QVAPORf48", "PRECIPf48_log10", "CLOUDf48_log10", "TCf48", "Pf48", "Uf48", "Vf48", "Wf48"])

for v in variables:
	data = np.fromfile("{}.bin.dat".format(v), dtype=np.float32)
	data = data.reshape([100, 500, 500])[:96, :496, :496]
	data.tofile("{}_truncated.bin.dat".format(v))

