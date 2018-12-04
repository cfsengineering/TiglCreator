//
// Created by cfse on 12/4/18.
//


#include "test.h" // Brings in the GTest framework
#include "tigl.h"

#include "CCPACSConfigurationManager.h"
#include "CCPACSConfiguration.h"
#include "CCPACSFuselage.h"

#include <string.h>
TEST(CreatorFuselage, getLength_SimpleModel)
{
const char* filename = "TestData/simpletest.cpacs.xml";

TiglCPACSConfigurationHandle tiglHandle = -1;
TixiDocumentHandle tixiHandle = -1;

ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
tigl::CCPACSConfiguration& config =  manager.GetConfiguration(tiglHandle);

tigl::CCPACSFuselage& fuselage =  config.GetFuselage(1);
double length = fuselage.GetLength();

ASSERT_NEAR(2, length, 0.0001);

tiglCloseCPACSConfiguration(tiglHandle);
tixiCloseDocument(tixiHandle);

}



TEST(CreatorFuselage, getLength_BoxWingModel)
{
    const char* filename = "TestData/boxWing.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config =  manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage =  config.GetFuselage(1);
    double length = fuselage.GetLength();

    ASSERT_NEAR(37, length, 2);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);

}


TEST(CreatorFuselage, getLength_CrmWingModel)
{
    const char* filename = "TestData/crm.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config =  manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage =  config.GetFuselage(1);
    double length = fuselage.GetLength();

    ASSERT_NEAR(61, length, 2);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);

}