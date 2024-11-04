"""
Function to read dataset from Lyn (1988)
"""


def readLynParams(pathParams, caseName=None):
    """
    pathData: str, path to folder where experimental data from Lyn are stored
    """
    caseList = []
    with open(pathParams, "r") as fData:
        for line in fData:
            if line[0] == "#":
                continue
            case = {}
            name, Hwater, Umean, uf, ws, dS = line[:-1].split(";")
            case["name"] = name
            case["Hwater"] = float(Hwater)
            case["Umean"] = float(Umean)
            case["uf"] = float(uf)
            case["ws"] = float(ws)
            case["dS"] = float(dS)
            if caseName == case["name"]:
                return case
            caseList.append(case)
    return caseList
