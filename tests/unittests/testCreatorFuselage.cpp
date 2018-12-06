//
// Created by cfse on 12/4/18.
//

#include "test.h" // Brings in the GTest framework
#include "tigl.h"

#include "CCPACSConfigurationManager.h"
#include "CCPACSConfiguration.h"
#include "CCPACSFuselages.h"
#include "CCPACSFuselage.h"

#include <string.h>

TEST(CreatorFuselage, getLength_SimpleModel)
{
    const char* filename = "TestData/simpletest.cpacs.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage = config.GetFuselage(1);
    double length                  = fuselage.GetLength();

    ASSERT_NEAR(2, length, 0.0001);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}

TEST(CreatorFuselage, getLength_BoxWingModel)
{
    const char* filename = "TestData/boxWing.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage = config.GetFuselage(1);
    double length                  = fuselage.GetLength();

    ASSERT_NEAR(37, length, 2);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}

TEST(CreatorFuselage, getLength_CrmWingModel)
{
    const char* filename = "TestData/crm.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage = config.GetFuselage(1);
    double length                  = fuselage.GetLength();

    ASSERT_NEAR(61, length, 2);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}


TEST(CreatorFuselage, getLengthBetween_SimpleModel)
{
    const char* filename = "TestData/simpletest.cpacs.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselage& fuselage = config.GetFuselage(1);
    double length                  = fuselage.GetLengthBetween("D150_Fuselage_1Section1IDElement1","D150_Fuselage_1Section2IDElement1");

    ASSERT_NEAR(1, length, 0.0001);

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}



TEST(CreatorFuselage, getLengthBetween_MultipleFuselagesModel)
{
    const char* filename = "TestData/multiple_fuselages.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);

    tigl::CCPACSFuselages& fuselages = config.GetFuselages();
    double r;
    std::string fuselageUID;

    fuselageUID = "SimpleFuselage";
    tigl::CCPACSFuselage& fuselage = fuselages.GetFuselage(fuselageUID);
    r = fuselage.GetLengthBetween("D150_Fuselage_1Section1IDElement1", "D150_Fuselage_1Section2IDElement1");
    EXPECT_DOUBLE_EQ(r, 1);

    r = fuselage.GetLengthBetween( "D150_Fuselage_1Section3IDElement1", "D150_Fuselage_1Section2IDElement1");
    EXPECT_DOUBLE_EQ(r, 1);


    fuselageUID = "SimpleFuselage4";
    tigl::CCPACSFuselage& fuselage2 = fuselages.GetFuselage(fuselageUID);
    r = fuselage2.GetLengthBetween("D150_Fuselage_4Section1IDElement1", "D150_Fuselage_4Section3IDElement1");
    EXPECT_DOUBLE_EQ(r, 4);

    // invalid input // TODO error management
    //EXPECT_THROW(fuselage2.GetLengthBetween("D150_Fuselage_4Section1IDElement1", "D150_Fuselage_4Section3IDElement1fas"), tigl::CTiglError  );

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}


TEST(CreatorFuselage, setLength_SimpleModel)
{
    const char* filename = "TestData/simpletest.cpacs.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config = manager.GetConfiguration(tiglHandle);

    double lengthBefore, lengthAfter;
    tigl::CCPACSFuselage& fuselage = config.GetFuselage(1);
    lengthBefore = fuselage.GetLength();
    EXPECT_NEAR(2,lengthBefore, 0.0001);
    fuselage.SetLength(5);
    lengthAfter = fuselage.GetLength();
    EXPECT_NEAR(5,lengthAfter, 0.0001);
    fuselage.SetLength(3.2);
    lengthAfter = fuselage.GetLength();
    EXPECT_NEAR(3.2, lengthAfter, 0.0001);

    ASSERT_EQ( SUCCESS, tixiSaveDocument(tixiHandle, "TestData/Output/simpletest-out.cpacs.xml" ) );

    lengthAfter = fuselage.GetLength();
    EXPECT_NEAR(3.2, lengthAfter, 0.0001);

    // TODO: write cpacs file from tigl
    /*
    fuselage.SetLengthBetween("D150_Fuselage_1Section1IDElement1","D150_Fuselage_1Section3IDElement1",10.5);
    lengthAfter = fuselage.GetLength();
    EXPECT_NEAR(10.5, lengthAfter, 0.0001);

    ASSERT_EQ( SUCCESS, tixiSaveDocument(tixiHandle, "TestData/Output/simpletest-out.cpacs.xml" ) );
    */

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);
}

TEST(CreatorFuselage, setLengthBetween_MultipleFuselagesModel){

    const char* filename = "TestData/multiple_fuselages.xml";

    TiglCPACSConfigurationHandle tiglHandle = -1;
    TixiDocumentHandle tixiHandle           = -1;

    ASSERT_EQ(SUCCESS, tixiOpenDocument(filename, &tixiHandle));
    ASSERT_EQ(TIGL_SUCCESS, tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

    tigl::CCPACSConfigurationManager& manager = tigl::CCPACSConfigurationManager::GetInstance();
    tigl::CCPACSConfiguration& config         = manager.GetConfiguration(tiglHandle);


    double newPartialL, oldPartialL,  globalL, oldGlobalL, r;
    std::string fuselageUID;

    tigl::CCPACSFuselages& fuselages = config.GetFuselages();

    fuselageUID = "SimpleFuselage";
    newPartialL = 3;
    tigl::CCPACSFuselage& fuselage1 = fuselages.GetFuselage(fuselageUID);
    fuselage1.SetLengthBetween("D150_Fuselage_1Section1IDElement1", "D150_Fuselage_1Section2IDElement1",  newPartialL);
    r = fuselage1.GetLengthBetween("D150_Fuselage_1Section1IDElement1", "D150_Fuselage_1Section2IDElement1");
    EXPECT_NEAR(r, newPartialL, 0.0001);

    // TODO need to implement rotation part of the setLengthBetween
    /*
    fuselageUID = "SimpleFuselage4";
    tigl::CCPACSFuselage& fuselage4 = fuselages.GetFuselage(fuselageUID);
    newPartialL = 5;
    fuselage4.SetLengthBetween("D150_Fuselage_4Section2IDElement1", "D150_Fuselage_4Section3IDElement1" , newPartialL);
    r = fuselage4.GetLengthBetween("D150_Fuselage_4Section2IDElement1", "D150_Fuselage_4Section3IDElement1");
    EXPECT_NEAR(r, newPartialL, 0.0001);
    */


    ASSERT_EQ( SUCCESS, tixiSaveDocument(tixiHandle, "TestData/Output/multiple-fuselages-out.cpacs.xml" ) );

    tiglCloseCPACSConfiguration(tiglHandle);
    tixiCloseDocument(tixiHandle);

}