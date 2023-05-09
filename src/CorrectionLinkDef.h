//
// Created by Misha on 3/7/2023.
//

#ifndef QNTOOLSINTERFACE_CORRECTIONLINKDEF_H
#define QNTOOLSINTERFACE_CORRECTIONLINKDEF_H

#if defined(__ROOTCLING__) || defined(__MAKECINT__)

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ class VariableManager;
#pragma link C++ class CorrectionTask;
#pragma link C++ class VectorConfig;
#pragma link C++ class TFilePtr;

#endif

#endif //QNTOOLSINTERFACE_CORRECTIONLINKDEF_H
