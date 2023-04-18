MPICH ERROR [Rank 1] [job id 7048372.0] [Thu Apr  6 22:00:06 2023] [nid004255] - Abort(2228111) (rank 1 in comm 0): Fatal error in PMPI_Test: Other MPI error, error stack:
PMPI_Test(190).................: MPI_Test(request=0x50638300, flag=0x7fff5923565c, status=0x7fff59235680) failed
MPIR_Test(75)..................: 
srun: error: nid004254: task 0: Segmentation fault
MPIR_Test_impl(36).............: 
MPIDI_Progress_test(89)........: 
MPIDI_OFI_handle_cq_error(1062): OFI poll failed (ofi_events.c:1064:MPIDI_OFI_handle_cq_error:Input/output error - VNI_NOT_FOUND)

aborting job:
Fatal error in PMPI_Test: Other MPI error, error stack:
PMPI_Test(190).................: MPI_Test(request=0x50638300, flag=0x7fff5923565c, status=0x7fff59235680) failed
MPIR_Test(75)..................: 
MPIR_Test_impl(36).............: 
MPIDI_Progress_test(89)........: 
MPIDI_OFI_handle_cq_error(1062): OFI poll failed (ofi_events.c:1064:MPIDI_OFI_handle_cq_error:Input/output error - VNI_NOT_FOUND)
srun: launch/slurm: _step_signal: Terminating StepId=7048372.0
srun: error: nid004255: task 1: Exited with exit code 255
srun: error: nid005335: task 47: Terminated
srun: error: nid005771: task 83: Terminated
srun: error: nid006209: task 105: Terminated
srun: error: nid006897: task 118: Terminated
srun: error: nid005834: task 92: Terminated
srun: error: nid006202: task 99: Terminated
srun: error: nid005546: task 52: Terminated
srun: error: nid006219: task 110: Terminated
srun: error: nid005702: task 73: Terminated
srun: error: nid006207: task 103: Terminated
srun: error: nid005636: task 61: Terminated
srun: error: nid005331: task 43: Terminated
srun: error: nid005334: task 46: Terminated
srun: error: nid006901: task 119: Terminated
srun: error: nid006205: task 101: Terminated
srun: error: nid004907: task 20: Terminated
srun: error: nid004661: task 16: Terminated
srun: error: nid005333: task 45: Terminated
srun: error: nid005767: task 79: Terminated
srun: error: nid005701: task 72: Terminated
srun: error: nid006997: task 122: Terminated
srun: error: nid005224: task 31: Terminated
srun: error: nid005553: task 59: Terminated
srun: error: nid005644: task 69: Terminated
srun: error: nid006201: task 98: Terminated
srun: error: nid004519: task 12: Terminated
srun: error: nid005547: task 53: Terminated
srun: error: nid005550: task 56: Terminated
srun: error: nid005833: task 91: Terminated
srun: error: nid007039: task 123: Terminated
srun: error: nid005828: task 86: Terminated
srun: error: nid005552: task 58: Terminated
srun: error: nid005223: task 30: Terminated
srun: error: nid005700: task 71: Terminated
srun: error: nid005328: task 40: Terminated
srun: error: nid005770: task 82: Terminated
srun: error: nid005336: task 48: Terminated
srun: error: nid005227: task 34: Terminated
srun: error: nid005226: task 33: Terminated
srun: error: nid005827: task 85: Terminated
srun: error: nid004911: task 24: Terminated
srun: error: nid005832: task 90: Terminated
srun: error: nid006204: task 100: Terminated
srun: error: nid005108: task 27: Terminated
srun: error: nid005330: task 42: Terminated
srun: error: nid006892: task 117: Terminated
srun: error: nid005261: task 39: Terminated
srun: error: nid006878: task 115: Terminated
srun: error: nid006212: task 107: Terminated
srun: error: nid007048: task 124: Terminated
srun: error: nid005637: task 62: Terminated
srun: error: nid004520: task 13: Terminated
srun: error: nid005643: task 68: Terminated
srun: error: nid006208: task 104: Terminated
srun: error: nid005337: task 49: Terminated
srun: error: nid005641: task 66: Terminated
srun: error: nid005765: task 77: Terminated
srun: error: nid005554: task 60: Terminated
srun: error: nid004664: task 17: Terminated
srun: error: nid004910: task 23: Terminated
srun: error: nid004908: task 21: Terminated
srun: error: nid005640: task 65: Terminated
srun: error: nid005230: task 37: Terminated
srun: error: nid005829: task 87: Terminated
srun: error: nid006218: task 109: Terminated
srun: error: nid007087: task 127: Terminated
srun: error: nid005518: task 51: Terminated
srun: error: nid006221: task 111: Terminated
srun: error: nid004623: task 15: Terminated
srun: error: nid005124: task 29: Terminated
srun: error: nid005123: task 28: Terminated
srun: error: nid005638: task 63: Terminated
srun: error: nid006242: task 112: Terminated
srun: error: nid005639: task 64: Terminated
srun: error: nid006245: task 113: Terminated
srun: error: nid004457: task 11: Terminated
srun: error: nid005329: task 41: Terminated
srun: error: nid005772: task 84: Terminated
srun: error: nid004913: task 26: Terminated
srun: error: nid004452: task 6: Terminated
srun: error: nid005902: task 94: Terminated
srun: error: nid005228: task 35: Terminated
srun: error: nid005835: task 93: Terminated
srun: error: nid004523: task 14: Terminated
srun: error: nid005768: task 80: Terminated
srun: error: nid004455: task 9: Terminated
srun: error: nid005231: task 38: Terminated
srun: error: nid004829: task 18: Terminated
srun: error: nid005225: task 32: Terminated
srun: error: nid005332: task 44: Terminated
srun: error: nid004906: task 19: Terminated
srun: error: nid007075: task 126: Terminated
srun: error: nid006752: task 114: Terminated
srun: error: nid005642: task 67: Terminated
srun: error: nid006884: task 116: Terminated
srun: error: nid004912: task 25: Terminated
srun: error: nid005548: task 54: Terminated
srun: error: nid006210: task 106: Terminated
srun: error: nid005766: task 78: Terminated
srun: error: nid004456: task 10: Terminated
srun: error: nid005549: task 55: Terminated
srun: error: nid004909: task 22: Terminated
srun: error: nid004453: task 7: Terminated
srun: error: nid004454: task 8: Terminated
srun: error: nid006968: task 121: Terminated
srun: error: nid006056: task 97: Terminated
srun: error: nid005764: task 76: Terminated
srun: error: nid005939: task 95: Terminated
srun: error: nid007049: task 125: Terminated
srun: error: nid005744: task 75: Terminated
srun: error: nid005830: task 88: Terminated
srun: error: nid005831: task 89: Terminated
srun: error: nid006911: task 120: Terminated
srun: error: nid006206: task 102: Terminated
srun: error: nid005769: task 81: Terminated
srun: error: nid005229: task 36: Terminated
srun: error: nid005645: task 70: Terminated
srun: error: nid005705: task 74: Terminated
srun: error: nid007088: task 128: Terminated
srun: error: nid005551: task 57: Terminated
srun: error: nid005517: task 50: Terminated
srun: error: nid004451: task 5: Terminated
srun: error: nid004450: task 4: Terminated
srun: error: nid004403: task 2: Terminated
srun: error: nid005940: task 96: Terminated
srun: error: nid004449: task 3: Terminated
srun: error: nid006216: task 108: Terminated
srun: Force Terminated StepId=7048372.0
