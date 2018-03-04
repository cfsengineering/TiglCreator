// Copyright (c) 2018 RISC Software GmbH
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

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <string>
#include <tixi.h>
#include "tigl_internal.h"

namespace tigl
{
class CTiglUIDManager;

namespace generated
{
    // This class is used in:
    // CPACSWindows

    // generated from /xsd:schema/xsd:complexType[245]
    class CPACSCutOut
    {
    public:
        TIGL_EXPORT CPACSCutOut(CTiglUIDManager* uidMgr);
        TIGL_EXPORT virtual ~CPACSCutOut();

        TIGL_EXPORT CTiglUIDManager& GetUIDManager();
        TIGL_EXPORT const CTiglUIDManager& GetUIDManager() const;

        TIGL_EXPORT virtual void ReadCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath);
        TIGL_EXPORT virtual void WriteCPACS(const TixiDocumentHandle& tixiHandle, const std::string& xpath) const;

        TIGL_EXPORT virtual const std::string& GetUID() const;
        TIGL_EXPORT virtual void SetUID(const std::string& value);

        TIGL_EXPORT virtual const boost::optional<std::string>& GetName() const;
        TIGL_EXPORT virtual void SetName(const boost::optional<std::string>& value);

        TIGL_EXPORT virtual const boost::optional<std::string>& GetDescription() const;
        TIGL_EXPORT virtual void SetDescription(const boost::optional<std::string>& value);

        TIGL_EXPORT virtual const double& GetWidth() const;
        TIGL_EXPORT virtual void SetWidth(const double& value);

        TIGL_EXPORT virtual const double& GetHeight() const;
        TIGL_EXPORT virtual void SetHeight(const double& value);

        TIGL_EXPORT virtual const double& GetFilletRadius() const;
        TIGL_EXPORT virtual void SetFilletRadius(const double& value);

        TIGL_EXPORT virtual const boost::optional<std::string>& GetReinforcementElementUID() const;
        TIGL_EXPORT virtual void SetReinforcementElementUID(const boost::optional<std::string>& value);

    protected:
        CTiglUIDManager* m_uidMgr;

        std::string                  m_uID;
        boost::optional<std::string> m_name;
        boost::optional<std::string> m_description;
        double                       m_width;
        double                       m_height;
        double                       m_filletRadius;
        boost::optional<std::string> m_reinforcementElementUID;

    private:
#ifdef HAVE_CPP11
        CPACSCutOut(const CPACSCutOut&) = delete;
        CPACSCutOut& operator=(const CPACSCutOut&) = delete;

        CPACSCutOut(CPACSCutOut&&) = delete;
        CPACSCutOut& operator=(CPACSCutOut&&) = delete;
#else
        CPACSCutOut(const CPACSCutOut&);
        CPACSCutOut& operator=(const CPACSCutOut&);
#endif
    };
} // namespace generated

// Aliases in tigl namespace
#ifdef HAVE_CPP11
using CCPACSCutOut = generated::CPACSCutOut;
#else
typedef generated::CPACSCutOut CCPACSCutOut;
#endif
} // namespace tigl
