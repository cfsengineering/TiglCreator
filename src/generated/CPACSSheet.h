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
#include <CCPACSPointXY.h>
#include <string>
#include <tixi.h>
#include "CPACSContinuityAtP.h"
#include "CreateIfNotExists.h"
#include "tigl_internal.h"

namespace tigl
{
class CTiglUIDManager;

namespace generated
{
    // This class is used in:
    // CPACSSheetList

    // generated from /xsd:schema/xsd:complexType[801]
    class CPACSSheet
    {
    public:
        TIGL_EXPORT CPACSSheet(CTiglUIDManager* uidMgr);
        TIGL_EXPORT virtual ~CPACSSheet();

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

        TIGL_EXPORT virtual const std::string& GetFromPointUID() const;
        TIGL_EXPORT virtual void SetFromPointUID(const std::string& value);

        TIGL_EXPORT virtual const boost::optional<CPACSContinuityAtP>& GetContinuityAtP1() const;
        TIGL_EXPORT virtual void SetContinuityAtP1(const boost::optional<CPACSContinuityAtP>& value);

        TIGL_EXPORT virtual const boost::optional<CCPACSPointXY>& GetOrientationAtP1() const;
        TIGL_EXPORT virtual boost::optional<CCPACSPointXY>& GetOrientationAtP1();

        TIGL_EXPORT virtual const std::string& GetToPointUID() const;
        TIGL_EXPORT virtual void SetToPointUID(const std::string& value);

        TIGL_EXPORT virtual const boost::optional<CPACSContinuityAtP>& GetContinuityAtP2() const;
        TIGL_EXPORT virtual void SetContinuityAtP2(const boost::optional<CPACSContinuityAtP>& value);

        TIGL_EXPORT virtual const boost::optional<CCPACSPointXY>& GetOrientationAtP2() const;
        TIGL_EXPORT virtual boost::optional<CCPACSPointXY>& GetOrientationAtP2();

        TIGL_EXPORT virtual CCPACSPointXY& GetOrientationAtP1(CreateIfNotExistsTag);
        TIGL_EXPORT virtual void RemoveOrientationAtP1();

        TIGL_EXPORT virtual CCPACSPointXY& GetOrientationAtP2(CreateIfNotExistsTag);
        TIGL_EXPORT virtual void RemoveOrientationAtP2();

    protected:
        CTiglUIDManager* m_uidMgr;

        std::string                         m_uID;
        boost::optional<std::string>        m_name;
        boost::optional<std::string>        m_description;
        std::string                         m_fromPointUID;
        boost::optional<CPACSContinuityAtP> m_continuityAtP1;
        boost::optional<CCPACSPointXY>      m_orientationAtP1;
        std::string                         m_toPointUID;
        boost::optional<CPACSContinuityAtP> m_continuityAtP2;
        boost::optional<CCPACSPointXY>      m_orientationAtP2;

    private:
#ifdef HAVE_CPP11
        CPACSSheet(const CPACSSheet&) = delete;
        CPACSSheet& operator=(const CPACSSheet&) = delete;

        CPACSSheet(CPACSSheet&&) = delete;
        CPACSSheet& operator=(CPACSSheet&&) = delete;
#else
        CPACSSheet(const CPACSSheet&);
        CPACSSheet& operator=(const CPACSSheet&);
#endif
    };
} // namespace generated

// Aliases in tigl namespace
#ifdef HAVE_CPP11
using CCPACSSheet = generated::CPACSSheet;
#else
typedef generated::CPACSSheet CCPACSSheet;
#endif
} // namespace tigl
