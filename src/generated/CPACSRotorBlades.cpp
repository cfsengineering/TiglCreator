// Copyright (c) 2016 RISC Software GmbH
//
// This file was generated by CPACSGen from CPACS XML Schema (c) German Aerospace Center (DLR/SC).
// Do not edit, all changes are lost when files are re-generated.
//
// Licensed under the Apache License, Version 2.0 (the "License")
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <CCPACSWing.h>
#include <cassert>
#include "CCPACSRotorcraftModel.h"
#include "TixiHelper.h"
#include "CTiglLogging.h"
#include "CTiglError.h"
#include "CPACSRotorBlades.h"

namespace tigl
{
    namespace generated
    {
        CPACSRotorBlades::CPACSRotorBlades(CCPACSRotorcraftModel* parent)
        {
            //assert(parent != NULL);
            m_parent = parent;
        }
        
        CPACSRotorBlades::~CPACSRotorBlades() {}
        
        CCPACSRotorcraftModel* CPACSRotorBlades::GetParent() const
        {
            return m_parent;
        }
        
        void CPACSRotorBlades::ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath)
        {
            // read element rotorBlade
            if (tixihelper::TixiCheckElement(tixiHandle, xpath + "/rotorBlade")) {
                tixihelper::TixiReadElements(tixiHandle, xpath + "/rotorBlade", m_rotorBlade, reinterpret_cast<CCPACSRotorBlades*>(this));
            }
            
        }
        
        void CPACSRotorBlades::WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const
        {
            // write element rotorBlade
            tixihelper::TixiSaveElements(tixiHandle, xpath + "/rotorBlade", m_rotorBlade);
            
        }
        
        const std::vector<unique_ptr<CCPACSWing> >& CPACSRotorBlades::GetRotorBlade() const
        {
            return m_rotorBlade;
        }
        
        std::vector<unique_ptr<CCPACSWing> >& CPACSRotorBlades::GetRotorBlade()
        {
            return m_rotorBlade;
        }
        
    }
}
