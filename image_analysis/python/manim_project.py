from manim import *

class PipelineOverview(Scene):
    def construct(self):
        # Title
        title = Text("Image Processing & Volume Estimation Pipeline", font_size = 36)
        title.to_edge(UP)
        self.play(Write(title))

        # Steps
        steps = [
            "1. Downscale & Segment Image (Ilastik)",
            "2. ROI Detection via Fiji",
            "3. Upscale & Clean Mask",
            "4. Skeletonize ROI into Midline",
            "5. Fit & Extend Midline Spline",
            "6. Compute Normals & Sample Widths",
            "7. Estimate Volume via Frustums",
            "8. Convert Volume to Mass"
        ]

        step_texts = VGroup(*[
            Text(step, font_size = 26).set_opacity(0.9)
            for step in steps
        ]).arrange(DOWN, aligned_edge = LEFT, buff = 0.4)

        step_texts.next_to(title, DOWN, buff = 0.6)

        # Animation
        self.play(LaggedStart(*[
            FadeIn(step, shift = RIGHT) for step in step_texts
        ], lag_ratio = 0.15))

        self.wait(2)

class IlastikSegmentation(Scene):
    def construct(self):
        
        original = ImageMobject("original.png").scale(0.75)
        downscaled = original.copy()
        segmentation = ImageMobject("segmented.png").scale(0.5).set_opacity(0.4)
   

        group = Group(original, downscaled).arrange(RIGHT, buff = 1)
        self.play(FadeIn(group))
        self.wait(1)
        
        segmentation.move_to(downscaled)
        label = Text("Segmentation via Ilastik", font_size = 24).next_to(downscaled, DOWN)

        self.play(FadeIn(segmentation), Write(label))
        self.wait(2)