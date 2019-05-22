import numpy as np

variables = np.array(["QC", "QG", "QI", "QR", "QS", "QV", "PRES", "RH", "T", "U", "V", "W"]) 

for v in variables:
	data = np.fromfile("{}-98x1200x1200.dat".format(v), dtype=np.float32)
	data = data.reshape([98, 1200, 1200])
	data = data[:96, :, :]
	if v.startswith('Q'):
		tmp = np.min(data[data > 0])
		data[data == 0] = tmp
		data = np.log10(data)
		data.tofile("{}_log10_truncated.bin.dat".format(v))
	else:
		data.tofile("{}_truncated.bin.dat".format(v))
