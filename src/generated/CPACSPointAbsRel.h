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

#pragma once

#include <tixi.h>
#include <string>
#include <boost/optional.hpp>
#include "tigl_internal.h"
#include <ECPACSTranslationType.h>

namespace tigl
{
    namespace generated
    {
        // This class is used in:
        // CPACSTransformation
        
        // generated from /xsd:schema/xsd:complexType[681]
        class CPACSPointAbsRel
        {
        public:
            TIGL_EXPORT CPACSPointAbsRel();
            TIGL_EXPORT virtual ~CPACSPointAbsRel();
            
            TIGL_EXPORT virtual void ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath);
            TIGL_EXPORT virtual void WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const;
            
            TIGL_EXPORT bool HasUID() const;
            TIGL_EXPORT const std::string& GetUID() const;
            TIGL_EXPORT void SetUID(const std::string& value);
            
            TIGL_EXPORT bool HasRefType() const;
            TIGL_EXPORT const ECPACSTranslationType& GetRefType() const;
            TIGL_EXPORT void SetRefType(const ECPACSTranslationType& value);
            
            TIGL_EXPORT bool HasX() const;
            TIGL_EXPORT const double& GetX() const;
            TIGL_EXPORT void SetX(const double& value);
            
            TIGL_EXPORT bool HasY() const;
            TIGL_EXPORT const double& GetY() const;
            TIGL_EXPORT void SetY(const double& value);
            
            TIGL_EXPORT bool HasZ() const;
            TIGL_EXPORT const double& GetZ() const;
            TIGL_EXPORT void SetZ(const double& value);
            
        protected:
            boost::optional<std::string>           m_uID;
            boost::optional<ECPACSTranslationType> m_refType;
            boost::optional<double>                m_x;
            boost::optional<double>                m_y;
            boost::optional<double>                m_z;
            
        private:
            #ifdef HAVE_CPP11
            CPACSPointAbsRel(const CPACSPointAbsRel&) = delete;
            CPACSPointAbsRel& operator=(const CPACSPointAbsRel&) = delete;
            
            CPACSPointAbsRel(CPACSPointAbsRel&&) = delete;
            CPACSPointAbsRel& operator=(CPACSPointAbsRel&&) = delete;
            #else
            CPACSPointAbsRel(const CPACSPointAbsRel&);
            CPACSPointAbsRel& operator=(const CPACSPointAbsRel&);
            #endif
        };
    }
    
    // This type is customized, use type CCPACSPointAbsRel
}
