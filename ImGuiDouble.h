#ifndef IMGUIDOUBLE_H
#define IMGUIDOUBLE_H

#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

// until libigl gets their act together and natively supports doubles
namespace ImGui
{

static bool DataTypeApplyOpFromTextScientific( const char* buf, const char* initial_value_buf, ImGuiDataType data_type, void* data_ptr, const char* scalar_format )
{
    scalar_format = "%f";
    float* v = (float*)data_ptr;
    const float old_v = *v;
    *v = (float)atof( buf );
    return *v != old_v;
}

static inline void DataTypeFormatString(ImGuiDataType data_type, void* data_ptr, const char* display_format, char* buf, int buf_size)
{
    if (data_type == ImGuiDataType_Int)
        ImFormatString(buf, buf_size, display_format, *(int*)data_ptr);
    else if (data_type == ImGuiDataType_Float)
        ImFormatString(buf, buf_size, display_format, *(float*)data_ptr);
}

bool InputScalarScientific( const char* label, ImGuiDataType data_type, void* data_ptr, const char* scalar_format, ImGuiInputTextFlags extra_flags )
{
    ImGuiWindow* window = GetCurrentWindow();
    if ( window->SkipItems )
        return false;

    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const ImVec2 label_size = CalcTextSize( label, NULL, true );

    ImGui::BeginGroup();
    ImGui::PushID( label );

    char buf[64];
    DataTypeFormatString( data_type, data_ptr, scalar_format, buf, IM_ARRAYSIZE( buf ) );

    bool value_changed = false;
    extra_flags |= ImGuiInputTextFlags_AutoSelectAll;
    if ( ImGui::InputText( "", buf, IM_ARRAYSIZE( buf ), extra_flags ) )
        value_changed = DataTypeApplyOpFromTextScientific( buf, GImGui->InputTextState.InitialText.begin(), data_type, data_ptr, scalar_format );

    ImGui::PopID();

    if ( label_size.x > 0 )
    {
        ImGui::SameLine( 0, style.ItemInnerSpacing.x );
        RenderText( ImVec2( window->DC.CursorPos.x, window->DC.CursorPos.y + style.FramePadding.y ), label );
        ItemSize( label_size, style.FramePadding.y );
    }
    ImGui::EndGroup();

    return value_changed;
}


bool InputFloatScientific( const char* label, float* v, const char *display_format = "%.3g", ImGuiInputTextFlags extra_flags = 0)
{
    return InputScalarScientific( label, ImGuiDataType_Float, (void*)v, display_format, extra_flags );
}

bool InputDoubleScientific( const char* label, double* v, const char *display_format = "%.3g", ImGuiInputTextFlags extra_flags = 0)
{
    float tmp = *v;
    bool ret = InputFloatScientific(label, &tmp, display_format, extra_flags);
    *v = tmp;
    return ret;
}

}

#endif
